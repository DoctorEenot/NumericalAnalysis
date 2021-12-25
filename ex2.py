import math
import random

def f(x):
    x_pow2 = pow(x,2)
    x_pow3 = x_pow2*x
    x_pow4 = x_pow3*x
    x_pow5 = x_pow4*x
    x_pow6 = x_pow5*x

    return 54*x_pow6 + 45*x_pow5 - 102*x_pow4\
        -69*x_pow3 + 35*x_pow2 + 16*x - 4

def df(x):
    x_pow2 = pow(x,2)
    x_pow3 = x_pow2*x
    x_pow4 = x_pow3*x
    x_pow5 = x_pow4*x

    return 324*x_pow5 + 225*x_pow4 - 408*x_pow3\
         - 207*x_pow2 + 70*x + 16
        
def d2f(x):
    x_pow2 = pow(x,2)
    x_pow3 = x_pow2*x
    x_pow4 = x_pow3*x

    return 1620*x_pow4 + 900*x_pow3 - 1224*x_pow2 - 414*x + 70



def Newton_Raphson(a,b,e,f,df,d2f):
    y1 = f(a)
    y2 = f(b)

    if y1*y2>=0:
        return (None,0)

    if(y1*d2f(a) > 0):
        x0 = a
    else:
        x0 = b


    try:
        xn = x0 - f(x0)/df(x0)
    except ZeroDivisionError:
        return (None,0)

    i = 0
    while math.fabs(x0-xn) > e:
        x0 = xn
        xn = x0 - (1/((df(x0)/f(x0)) - 0.5*(d2f(x0)/df(x0))))
        i += 1
    
    return (round(xn,5),i)


def dichotomy_method(a,b,e,f):
    y1 = f(a)
    
    y3 = 1
    i = 0
    while math.fabs(y3) > e:
        x = random.uniform(a,b)
        y3 = f(x)
        if y1*y3 < 0:
            # if y1 and y3 have different signes
            b = x
        else:
            a = x
        i += 1

    return (round(x,5),i)


def secant_method(xn,xn_1,xn_2,e,f):
    def q():
        return f(xn)/f(xn_1)
    def r():
        return f(xn_2)/f(xn_1)
    def s():
        return f(xn_2)/f(xn)

    # xn_3 = xn_2 - ((r()*(r()-q()) * (xn_2-xn_1)) + ((1-r())*s()*(xn_2-xn)))/((q()-1)*(r()-1)*(s()-1))
    # y = f(xn_3)
    
    # xn_2 = xn_1
    # xn_1 = xn
    # xn = xn_3
    y = 1

    i = 0
    while math.fabs(y) > e:
        xn_3 = xn_2 - ((r()*(r()-q()) * (xn_2-xn_1)) + ((1-r())*s()*(xn_2-xn)))/((q()-1)*(r()-1)*(s()-1))
        y = f(xn_3)

        xn_2 = xn_1
        xn_1 = xn
        xn = xn_3
        i += 1

    return (round(xn,5),i)

        


def main():

    result = Newton_Raphson(-2,-1,0.00001,f,df,d2f)
    print("X1 with Newton-Raphson:",result[0],"Repetitions:",result[1])

    result = Newton_Raphson(-1,2,0.00001,f,df,d2f)
    print("X2 with Newton-Raphson:",result[0],"Repetitions:",result[1])

    result = Newton_Raphson(0.0,0.3,0.00001,f,df,d2f)
    print("X3 with Newton-Raphson:",result[0],"Repetitions:",result[1])

    result = Newton_Raphson(0.4,0.6,0.00001,f,df,d2f)
    print("X4 with Newton-Raphson:",result[0],"Repetitions:",result[1])

    result = Newton_Raphson(1,2,0.00001,f,df,d2f)
    print("X5 with Newton-Raphson:",result[0],"Repetitions:",result[1])

    print()

    result = dichotomy_method(-2,-1,0.00001,f)
    print("X1 with Dichotomy:",result[0],"Repetitions:",result[1])

    result = dichotomy_method(-1,-0.66667,0.00001,f)
    print("X2 with Dichotomy:",result[0],"Repetitions:",result[1])

    result = dichotomy_method(0,0.3,0.00001,f)
    print("X3 with Dichotomy:",result[0],"Repetitions:",result[1])

    result = dichotomy_method(0.4,0.6,0.00001,f)
    print("X4 with Dichotomy:",result[0],"Repetitions:",result[1])

    result = dichotomy_method(1,2,0.00001,f)
    print("X5 with Dichotomy:",result[0],"Repetitions:",result[1])

    print()

    result = secant_method(-4,-3,-2,0.00001,f)
    print("X1 with Secant:",result[0],"Repetitions:",result[1])

    result = secant_method(-3,-2,-1,0.00001,f)
    print("X2 with Secant:",result[0],"Repetitions:",result[1])

    result = secant_method(0,0.2,0.3,0.00001,f)
    print("X3 with Secant:",result[0],"Repetitions:",result[1])

    result = secant_method(0.3,0.4,0.5,0.00001,f)
    print("X4 with Secant:",result[0],"Repetitions:",result[1])

    result = secant_method(1,2,3,0.00001,f)
    print("X5 with Secant:",result[0],"Repetitions:",result[1])


if __name__ == '__main__':
    main()