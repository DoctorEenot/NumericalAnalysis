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

def create_least_squares(x,y):
    if len(x) != len(y):
        raise Exception("Array lengths are different")

    coefitients = []
    for i in range(len(x)-1):
        x0 = x[i]
        x1 = x[i+1]
        y0 = y[i]
        y1 = y[i+1]
        x_sum = x0 + x1
        y_sum = y0 + y1
        xy_sum = (x0*y0) + (x1*y1)
        x_square_sum = (x0*x0) + (x1*x1)
        
        a = ((2*xy_sum) - (x_sum*y_sum))\
                /((2*x_square_sum) - (x_sum*x_sum))

        b = (y_sum - (a*x_sum))/2

        coefitients.append((a,b))

    def least_squares(x_new):
        index = min(bisect.bisect(x, x_new)-1, len(x)-2)
        a,b = coefitients[index]
        y = (a*x_new) + b
        return y

    return least_squares
        


def main():
    TESTING_BOUNDARIES = (-pi,pi)
    x_step = 0.1
    x = TESTING_BOUNDARIES[0]
    digits = 200
    x_values = []
    y_lib = []
    y_approximated = []
    print("**TESTING SIN APPROXIMATION WITH POLYNOMIALS**\n")

    error_sum = 0
    count = 0
    while x<=TESTING_BOUNDARIES[1]:
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
    while x<=TESTING_BOUNDARIES[1]:
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
    least_squares = create_least_squares(x_values,y_lib)
    x_step = 0.001
    x_values = []

    error_sum = 0
    count = 0
    while x<=TESTING_BOUNDARIES[1]:
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
