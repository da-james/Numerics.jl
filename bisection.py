import math

def f(x):
    return math.sqrt(x) - math.cos(x)

def bisection(start, end):
    err = .0001
    a = start
    b = end
    p = (a+b)/2
    print("a:",a)
    print("b:",b)
    print("p:",p)
    while(f(p) < -1*err or f(p) > err):
        print("f(a) =", f(a))
        print("f(b) =", f(b))
        print("f(p) =", f(p))
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
    else:
        print("found p is",p)


bisection(0,1)
