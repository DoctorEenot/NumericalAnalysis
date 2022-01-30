import numpy as np

pow = np.power

def create_polynomial(x_values,y_values,polynomial_degree=2):
    '''
        Creates matrixes to be solved
    '''

    if polynomial_degree<1:
        raise Exception("Bad degree")

    coefs = np.zeros((polynomial_degree+1,polynomial_degree+1),
                        dtype=np.float)

    results = np.zeros((polynomial_degree+1,),dtype=np.float)


    c_sum = len(x_values)
    coefs[polynomial_degree][polynomial_degree] = c_sum
    for i in range(1,polynomial_degree+1):
        sm = sum(pow(x,i) for x in x_values)
        coefs[polynomial_degree][polynomial_degree-i] = sm

    y_sum = sum(y_values)
    results[polynomial_degree] = y_sum

    c = 1
    for count in range(polynomial_degree-1,-1,-1):
        for i in range(0,polynomial_degree+1):
            sm = sum(pow(x,c+i) for x in x_values)
            coefs[count][polynomial_degree-i] = sm

        y_sum = sum(y*pow(x,c) for x,y in zip(x_values,y_values))
        results[count] = y_sum  

        c += 1          


    return (coefs,results)

def get_least_square_solver(polynomial,result):
    degree = polynomial.shape[0]
    coefficients = np.flip(np.linalg.solve(polynomial,result))

    def approximator(x:float):
        y = coefficients[0]
        for k in range(1,degree):
            y += coefficients[k]*pow(x,k)

        return y

    return approximator




def main():
    x_values = np.array([2,3,4,5],dtype=np.float)
    y_values = np.array([2,2,2,2],dtype=np.float)
    polynomial,res = create_polynomial(x_values,y_values,2)
    print(polynomial)
    print(res)

    solver = get_least_square_solver(polynomial,res)
    print(solver(6.0))





if __name__ == "__main__":
    main()