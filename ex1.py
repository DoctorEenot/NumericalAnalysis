import matplotlib.pyplot as plt 
import math


def f(x):
    e_pow = pow(math.e,x-2)
    x_pow_2 = pow(x,2)
    x_pow_3 = x_pow_2*x
    return 14*x*e_pow-12*e_pow-7*x_pow_3+20*x_pow_2-26*x+12

# f'
def df(x):
    e_pow = pow(math.e,x-2)
    return 14*x*e_pow+2*e_pow-21*pow(x,2)+40*x-26

# f''
def d2f(x):
    e_pow = pow(math.e,x-2)
    return 14*x*e_pow+16*e_pow-42*x+40


def dichotomy_method(a,b,e,f):
    y1 = f(a)
    y2 = f(b)

    if y1*y2>=0:
        return (None,0)
    
    x = (a+b)/2
    y3 = f(x)
    i = 0
    while math.fabs(y3) > e:
        x = (a+b)/2
        y3 = f(x)
        if y1 * y3 < 0:
            # if y1 and y3 have different signes
            b = x
        else:
            a = x
        i += 1

    return (round(x,5),i)

def Newton_Raphson_method(a,b,e,f,df,d2f):
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
        xn = x0 - f(x0)/df(x0)
        i += 1
    
    return (round(xn,5),i)


def secand_method(x,e,f,df):
    #x1 = x0 + 0.1

    y = f(x)
    i = 0
    
    while math.fabs(y) > e:
        h = -(y/df(x))
        x += h
        y = f(x)
        i += 1

    return (round(x,5),i)




def print_graph(a,b,step):
    xpoints = []
    ypoints = []

    x = a
    x_step = step

    while x <= b:
        xpoints.append(x)
        ypoints.append(f(x))

        x += x_step

    plt.plot(xpoints, ypoints)
    plt.show()





def main():
    print_graph(0,3,0.1)

    result = dichotomy_method(0,1,0.00001,f)
    print("X1 with dichotomy:",result[0],"Repetitions:",result[1])
    result = dichotomy_method(1,10,0.00001,f)
    print("X2 with dichotomy:",result[0],"Repetitions:",result[1])

    print()

    result = Newton_Raphson_method(0,1,0.00001,f,df,d2f)
    print("X1 with Newton-Raphson:",result[0],"Repetitions:",result[1])
    result = Newton_Raphson_method(1,10,0.00001,f,df,d2f)
    print("X2 with Newton-Raphson:",result[0],"Repetitions:",result[1])

    print()

    result = secand_method(0,0.00001,f,df)
    print("X1 with Newton-Raphson:",result[0],"Repetitions:",result[1])
    result = secand_method(10,0.00001,f,df)
    print("X2 with Newton-Raphson:",result[0],"Repetitions:",result[1])


if __name__ == "__main__":
    main()


