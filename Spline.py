import bisect

def changes(x):
    return [x[i+1] - x[i] for i in range(len(x) - 1)]

def create_tridiagonalmatrix(n: int, h):
    A = [h[i] / (h[i] + h[i + 1]) for i in range(n - 2)] + [0]
    B = [2] * n
    C = [0] + [h[i + 1] / (h[i] + h[i + 1]) for i in range(n - 2)]
    return A, B, C

def create_target(n: int, h, y):
    return [0] + [6 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]) / (h[i] + h[i-1]) for i in range(1, n - 1)] + [0]

def solve_tridiagonalsystem(A, B, C, D):
    c_p = C + [0]
    d_p = [0] * len(B)
    X = [0] * len(B)

    c_p[0] = C[0] / B[0]
    d_p[0] = D[0] / B[0]
    for i in range(1, len(B)):
        c_p[i] = c_p[i] / (B[i] - c_p[i - 1] * A[i - 1])
        d_p[i] = (D[i] - d_p[i - 1] * A[i - 1]) / (B[i] - c_p[i - 1] * A[i - 1])

    X[-1] = d_p[-1]
    for i in range(len(B) - 2, -1, -1):
        X[i] = d_p[i] - c_p[i] * X[i + 1]

    return X

def create_spline(x, y):
    length_x = len(x)

    if length_x != len(y):
        raise Exception('Array lengths are different')

    computed_changes = changes(x)

    A, B, C = create_tridiagonalmatrix(length_x, computed_changes)
    D = create_target(length_x, computed_changes, y)

    M = solve_tridiagonalsystem(A, B, C, D)

    coefficients = [[(M[i+1]-M[i])*computed_changes[i]*computed_changes[i]/6, 
                    M[i]*computed_changes[i]*computed_changes[i]/2, 
                    (y[i+1] - y[i] - (M[i+1]+2*M[i])*computed_changes[i]*computed_changes[i]/6), 
                    y[i]] for i in range(length_x-1)]

    def spline(val):
        idx = min(bisect.bisect(x, val)-1, length_x-2)
        z = (val - x[idx]) / computed_changes[idx]
        C = coefficients[idx]
        return (((C[0] * z) + C[1]) * z + C[2]) * z + C[3]

    return spline