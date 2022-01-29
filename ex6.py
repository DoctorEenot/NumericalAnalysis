from math import sin,pi



def trapezodial(x_values:list,f):
    result = 0.0

    if len(x_values)<2:
        raise Exception("Not enough x values")


    for i in range(0,len(x_values)-1):
        a = x_values[i]
        b = x_values[i+1]

        result += ((b-a)/2)*(f(a)+f(b))


    return result

def Simpson(x_values,f):
    result = 0.0

    if len(x_values)<2:
        raise Exception("Not enough x values")

    a = x_values[0]
    b = x_values[1]

    result += ((b-a)/6)*(f(a)+f(b))

    for i in range(1,len(x_values)-2):
        a = x_values[i]
        b = x_values[i+1]
        
        if i%2 == 0:
            result += ((b-a)/6)*(f(a)+f(b))*4
        else:
            result += ((b-a)/6)*(f(a)+f(b))*2

    if len(x_values)>2:
        a = x_values[len(x_values)-2]
        b = x_values[len(x_values)-1]

        result += ((b-a)/6)*(f(a)+f(b))


    return result

def percent_error(e_value,r_value):
    return (abs(e_value-r_value)/r_value)*100

def print_results(result,real_result):
    print("- Approximated result:",result)

    p_error = percent_error(result,real_result)
    print(f"- Percent error: {p_error}%")

    numerical_error = real_result-result
    print("- Numerical error:",numerical_error)

def main():
    N = 11
    TESTING_BOUNDARIES = (0,pi/2)
    REAL_RESULT = 1

    # Preparing testing data
    x_values = []
    x_value = TESTING_BOUNDARIES[0]
    step = (TESTING_BOUNDARIES[1]-TESTING_BOUNDARIES[0])/N

    while x_value<TESTING_BOUNDARIES[1]-step:
        x_values.append(x_value)
        x_value += step

    
    x_values.append(TESTING_BOUNDARIES[1])

    print("**TESTING TRAPIZODIAL RULE**\n")

    result = trapezodial(x_values,sin)
    print_results(result,REAL_RESULT)

    print()
    print("**TESTING SIMPSON'S RULE**\n")

    result = Simpson(x_values,sin)
    
    print_results(result,REAL_RESULT)





if __name__ == "__main__":
    main()    
        

