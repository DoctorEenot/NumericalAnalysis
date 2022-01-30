from math import factorial,sin,pi
import numpy as np
from decimal import *
import matplotlib.pyplot as plt
import Spline
import bisect


def sin_polynomial(x:float,digits:int):
    y = Decimal(x)
    getcontext().prec = digits
    step = 0
    n = 3
    while step < digits:
        y += Decimal(pow(-1,n) * pow(x,n)) / Decimal(factorial(n))
        n += 2
        step += 1

    return y

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
    TESTING_BOUNDARIES = (-pi,pi)
    x_step = 0.01
    x = TESTING_BOUNDARIES[0]
    digits = 200
    x_values = []
    y_lib = []
    y_approximated = []
    print("**TESTING SIN APPROXIMATION WITH POLYNOMIALS**\n")

    error_sum = 0
    count = 0
    while x<=TESTING_BOUNDARIES[1] and\
            count<digits:
        x_values.append(x)
        sin_x = Decimal(sin(x))
        y_lib.append(sin_x)

        #print(f"- From Math sin({x}) =",sin_x)

        sin_x_approximated = sin_polynomial(x,digits)
        y_approximated.append(sin_x_approximated)
        #print(f"- Approximated ({digits} digits) sin({x}) =",sin_x_approximated)
        #print("- Error: ",abs(sin_x-sin_x_approximated))
        error_sum += abs(sin_x-sin_x_approximated)
        count += 1
        x += x_step

    avg_error = error_sum/count
    print("- Average Error:",avg_error)
    
    plt.plot(x_values, y_lib, color='r', label='sin lib')
    plt.plot(x_values, y_approximated, color='g', label='sin aproximate')
    plt.legend()
    plt.show()

    print()
    print("**TESTING SIN APPROXIMATION WITH SPLINES**\n")

    x_values = []
    y_lib = []
    x_step = 1
    y_approximated = []
    x = TESTING_BOUNDARIES[0]
    while x<=TESTING_BOUNDARIES[1]:
        x_values.append(x)
        sin_x = sin(x)
        y_lib.append(sin_x)
        x += x_step
    
    plt.plot(x_values, y_lib, color='r', label='sin lib')

    x = TESTING_BOUNDARIES[0]
    spline = Spline.create_spline(x_values,y_lib)
    x_step = 0.01
    x_values = []

    error_sum = 0
    count = 0
    while x<=TESTING_BOUNDARIES[1] and\
            count<digits:
        sin_x_approximated = spline(x)
        y_approximated.append(sin_x_approximated)
        x_values.append(x)

        sin_x = sin(x)
        error_sum += abs(sin_x-sin_x_approximated)
        x += x_step
        count += 1

    avg_error = error_sum/count
    print("- Average Error:",avg_error)
    
    plt.plot(x_values, y_approximated, color='g', label='sin aproximate')
    plt.legend()
    plt.show()

    print()
    print("**TESTING SIN APPROXIMATION WITH LEAST SQUARES**\n")

    x_values = []
    y_lib = []
    x_step = 0.1
    y_approximated = []
    x = TESTING_BOUNDARIES[0]

    while x<=TESTING_BOUNDARIES[1]:
        x_values.append(x)
        sin_x = sin(x)
        y_lib.append(sin_x)
        x += x_step
    
    plt.plot(x_values, y_lib, color='r', label='sin lib')

    x = TESTING_BOUNDARIES[0]
    polynomial,res = create_polynomial(x_values,y_lib,3)
    least_squares = get_least_square_solver(polynomial,res)
    x_step = 0.01
    x_values = []

    error_sum = 0
    count = 0
    while x<=TESTING_BOUNDARIES[1] and\
            count<digits:
        sin_x_approximated = least_squares(x)
        y_approximated.append(sin_x_approximated)
        x_values.append(x)
        
        sin_x = sin(x)
        error_sum += abs(sin_x-sin_x_approximated)
        x += x_step
        count += 1

    avg_error = error_sum/count
    print("- Average Error:",avg_error)
    
    plt.plot(x_values, y_approximated, color='g', label='sin aproximate')
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
