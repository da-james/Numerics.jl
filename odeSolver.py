from numpy import zeros,log,ones,diag,linalg,arange

def aitkenEst(yh,yh2,yh4):
    return log((yh-yh2)/(yh2-yh4))/log(2)

def solveODE(f,ti,tf,yi,funcName,n=0,h=0):
    if(h == 0 and n != 0):
        h = (tf-ti)/n
    elif(h != 0 and n == 0):
        n = int((tf-ti)/h)
        
    u = zeros((2,n+1))
    
    u[0,0] = ti
    u[1,0] = yi
    
    for i in range(n):
        u[0,i+1] = u[0,i] + h
        if(funcName == 'euler'):
            u[1,i+1] = euler(f,u[0,i],u[1,i],h)
        elif(funcName == 'rk2'):
            u[1,i+1] = rk2(f,u[0,i],u[1,i],h)
        elif(funcName == 'rk4'):
            u[1,i+1] = rk4(f,u[0,i],u[1,i],h)
        elif(funcName == 'trapezoid'):
            u[1,i+1] = trapezoid(f,u[0,i],u[1,i],h)
        elif(funcName == 'ab2'):
            if(i == 0): u[1,i+1] = rk4(f,u[0,i],u[1,i],h)
            else: u[1,i+1] = ab2(f,u[0,i],u[1,i],h,u[0,i-1],u[1,i-1])
        elif(funcName == 'am2'):
            if(i == 0): u[1,i+1] = rk4(f,u[0,i],u[1,i],h)
            else: u[1,i+1] = am2(f,u[0,i],u[1,i],h,u[0,i-1],u[1,i-1])
        else:
            print('no method was chosen')
            break
            
    return u

def solveSYS(f,ti,tf,yi,m, n=0, h=0):
    if(h == 0 and n != 0):
        h = (tf-ti)/n
    elif(h != 0 and n == 0):
        n = int((tf-ti)/h)
        
    t = zeros((n+1))
    y = zeros((m,n+1))
    
    t[0] = ti
    y[:,0] = yi
    
    for i in range(n):
        t[i+1] = t[i] + h
        
        
        s1 = h * f(t[i],y[:,i])
        s2 = h * f(t[i] + h/2, y[:,i] + s1/2)
        
        y[:,i+1] = y[:,i] + s2
        
    return y,t   
    
    
def solveODE2(f,p,q,ta,ya,tb,yb,h=0,n=0):
    if(h == 0 and n != 0):
        h = (tb-ta)/n
    elif(h != 0 and n == 0):
        n = int((tb-ta)/h)
        
    m = n-1
    
    s1 = 1/h**2 - p/(2*h)
    s2 = q - 2/h**2
    s3 = 1/h**2 + p/(2*h)
    
    x = arange(ta,tb+h,h)
        
    L = diag(ones(m-1),-1) * s1
    D = diag(ones(m))      * s2
    U = diag(ones(m-1),1)  * s3
    A = L + D + U
        
    f_b = f(x[1:n])
    f_b[0]  -= s1*ya
    f_b[-1] -= s3*yb
    
    y = ones(x.size)
    y[0] = ya
    y[-1] = yb
    y[1:n] = linalg.solve(A,f_b)
    
    return y,x   
    
    
def euler(f,ti,yi,h):
    y = yi + h*f(ti, yi)
    
    return y

def rk2(f,ti,yi,h):
    s1 = f(ti, yi)
    s2 = f(ti + h, yi + h*s1)
    y = yi + h/2*(s1 + s2)
        
    return y

def rk4(f,ti,yi,h):
    s1 = h*f(ti, yi)
    s2 = h*f(ti + h/2, yi + s1/2)
    s3 = h*f(ti + h/2, yi + s2/2)
    s4 = h*f(ti + h, yi + s3)
    
    y = yi + (s1 + 2*s2 + 2*s3 + s4)/6
        
    return y

def trapezoid(f,ti,yi,h):
    s1 = rk4(f,ti,yi,h)
    y = yi + h/2*(f(ti,yi) + s1)
    
    return y

def ab2(f,ti,yi,h,tf,yf):
    y = yi + 3/2*h*f(ti,yi) - 1/2*h*f(tf,yf)
    
    return y

def am2(f,ti,yi,h,tf,yf):
    s1 = rk4(f,ti,yi,h)
    
    y = yi + h/12*(5*s1+8*f(ti,yi)-f(tf,yf))
    
    return y
    