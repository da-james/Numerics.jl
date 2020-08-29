import math

def f(x):
    return math.sqrt(x) - math.cos(math.radians(x))
    #return x**2 - 3

def bisection(error, start, end):
    err = error
    n = 0
    a = start
    b = end
    p = (a+b)/2
    shift = 0
    leftShift = shift-1*err
    rightShift = shift+err
    
    print("a:",a)
    print("b:",b)
    print("p:",p)
    
    while(f(p) < leftShift or f(p) > rightShift):
        print("f(a) =", f(a))
        print("f(b) =", f(b))
        print("f(p) =", f(p))
        n += 1
        if(f(p) > 0 and f(a) < 0 ):
            b = p
            print("new b:",b)
        elif(f(p) > 0 and f(b) < 0):
            a = p
            print("new a:",a)
        elif(f(p) < 0 and f(a) > 0):
            b = p
            print("new b:",b)
        elif(f(p) < 0 and f(b) > 0):
            a = p
            print("new a:",a)
        else:
            print("exception made")
            break
        p = (a+b)/2
        print("This is step:",n)
        print("p =", p)
    else:
        print("found p is",p)
        print("Steps:",n)


print("Bisection method will try to find a zero of a given function.")
print("Please input the following: ")
error = float(input("Error: "))
leftB = int(input("Left Bound: "))
rightB = int(input("Right Bound: "))

bisection(error,leftB,rightB)
