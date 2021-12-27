import numpy as np
from math import sqrt


def lu_decomposition(A):
    n = A.shape[0]
    
    U = A.copy()
    L = np.eye(n, dtype=np.double)
    
    for i in range(n):
        factor = U[i+1:, i] / U[i, i]
        L[i+1:, i] = factor
        U[i+1:] -= factor[:, np.newaxis] * U[i]
        
    return L, U

def forward_substitution(L, b):
    n = L.shape[0]
    
    y = np.zeros_like(b, dtype=np.double)
    
    y[0] = b[0] / L[0, 0]

    for i in range(1, n):
        dotted = np.dot(L[i,:i], y[:i])
        y[i] = (b[i] - dotted) / L[i,i]
        
    return y

def back_substitution(U, y):
    n = U.shape[0]
    
    x = np.zeros_like(y, dtype=np.double)

    x[-1] = y[-1] / U[-1, -1]

    for i in range(n-2, -1, -1):
        dotted = np.dot(U[i,i:], x[i:])
        x[i] = (y[i] - dotted) / U[i,i]
        
    return x



def plu_solve(A,b):
    L, U = lu_decomposition(A)
    y = forward_substitution(L, b)   
    x = back_substitution(U, y)
    return x

def cholesky(A):
    n = A.shape[0]

    L = np.zeros((n,n),dtype=np.double)

    for i in range(n):
        for k in range(i+1):
            tmp_sum = sum(L[i][j] * L[k][j] for j in range(k))
            
            if (i == k): 
                L[i][k] = sqrt(A[i][i] - tmp_sum)
            else:
                L[i][k] = (1.0 / L[k][k] * (A[i][k] - tmp_sum))
    return L


def gen_matrix(n):
    '''
        Generates matrix with 
        A(i; i) = 5, A(i + 1; i) = A(i; i + 1) = -2
    '''
    A = np.zeros(shape=(n,n),dtype=np.double)

    A[0][0] = 5
    A[0][1] = -2
    A[n-1][n-1] = 5
    A[n-1][n-2] = -2

    for i in range(1,n-1):
        A[i][i] = 5
        A[i][i-1] = -2
        A[i][i+1] = -2
    
    return A

def GaussSeidel(A, b, eps):
    n = A.shape[0]
    x = np.zeros(shape=n,dtype=np.double) 

    found = False
    while not found:
        x_buf = np.copy(x)
        for i in range(n):
            sum_1 = sum(A[i][j] * x_buf[j] for j in range(i))
            sum_2 = sum(A[i][j] * x[j] for j in range(i + 1, n))
            x_buf[i] = (b[i] - sum_1 - sum_2) / A[i][i]

        found = np.sqrt(sum((x_buf[i] - x[i]) ** 2 for i in range(n))) <= eps
        x = x_buf

    return x    

def main():
    A = np.array([[7, 3, -1, 2], 
                [3, 8, 1, -4], 
                [-1, 1, 4, -1], 
                [2, -4, -1, 6]],dtype=np.double)
    b = np.array([[1,1,1,1]],dtype=np.double).T
    
    x = plu_solve(A,b)
    print("Result with PLU method:")
    print(x)

    print()

    L = cholesky(A)
    print("Cholesky L:")
    print(L)

    print()

    print("Gauss-Seidel on matrix with n=10")
    A = gen_matrix(10)

    # generating b matrix
    b = np.zeros(shape=(10,1),dtype=np.double)
    b[0][0] = 3
    b[9][0] = 3
    for i in range(1,9):
        b[i][0] = 1

    x = GaussSeidel(A,b,0.0001)
    print("Result:")
    print(x)

    print()

    print("Gauss-Seidel on matrix with n=10000")
    A = gen_matrix(10000)

    # generating b matrix
    b = np.zeros(shape=(10000,1),dtype=np.double)
    b[0][0] = 3
    b[10000-1][0] = 3
    for i in range(1,10000-1):
        b[i][0] = 1

    x = GaussSeidel(A,b,0.0001)
    print("Result:")
    print(x)


    
    



if __name__ == "__main__":
    main()