import math as ma

def f(x):
    return 4*x**3 + 2*x - 2;

def fprime(x):
    return 12*x**2 + 2;

def newton(error,left):
    # Intial Values
    n = 0;
    p = left;
    nextP = p;
    cond = True;
    absErr = ma.fabs(p - nextP);

    while(cond):
        print("Step:",n);
        print("p =", p);
        print("Error: ",absErr);

        nextP = p - (f(p)/fprime(p));
        absErr = ma.fabs(p-nextP);
        p = nextP;
        n += 1;

        if(absErr < error):
            cond = False;

    # final print for the last values
    print("Step:",n);
    print("p =",p);
    print("Error:",absErr);

def secant(error,left,right,pZero,pOne):
    # Intial Values
    n = 1;
    prevP = pZero;
    p = pOne;
    nextP = p;
    
    cond = True;
    
    absErr = ma.fabs(p - nextP);

    while(cond):
        print("Step:",n);
        print("p =", p);
        print("Error: ",absErr);

        prevQ = f(prevP);
        q = f(p);

        nextP = p - q*((p-prevP)/(q-prevQ));
        absErr = ma.fabs(p-nextP);
        prevP = p;
        p = nextP;
        n += 1;

        if(absErr < error):
            cond = False;

    # final print for the last values
    print("Step:",n);
    print("p =",p);
    print("Error:",absErr);


print("The Newton Iteration will try to find a zero of a function.");
error = .0001;
leftB = 1;
# rightB = 1;

newton(error,leftB);
# print("****************************************")
# print("The Secant Iteration will now also try to find a zero of the function.")
# secant(error,leftB,rightB,1,1.5)
