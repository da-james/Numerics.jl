import math as ma

def g(x):
    return .5*(3/x + x);

def fixedPoint(error, start, end):
    # Intial Values
    # step
    n = 0;
    # starting p
    p = float((start+end)/2);
    # conditional for while loop
    cond = True;
    # next p value to be calculated
    nextP = p;
    # error between p values
    absErr = float(ma.fabs(p - nextP));
    
    while(cond):
        # prints out current step and it's value
        print("Step:",n);
        print("p =", p);
        print("Error: ",absErr);

        # calculates next p
        nextP = g(p);
        # calculates error
        absErr = float(ma.fabs(p - nextP));
        # sets the new p to the current p
        p = nextP;
        # increments the step
        n += 1;

        # if the error is less than the given error
        # then break out of the loop
        if(absErr < error):
            cond = False;

    # final print for the last values
    print("Step:",n);
    print("p =",p);
    print("Error:",absErr);


            
print("The Fixed Point Iteration will try to find a zero of a function.");
error = .0001;
leftB = 1;
rightB = 2;

fixedPoint(error,leftB,rightB);
