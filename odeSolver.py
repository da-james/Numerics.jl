from numpy import zeros,log

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
        elif(funcName == 'heun'):
            u[1,i+1] = heun(f,u[0,i],u[1,i],h)
        elif(funcName == 'rk4'):
            u[1,i+1] = rk4(f,u[0,i],u[1,i],h)
        elif(funcName == 'trapezoid'):
            u[1,i+1] = trapezoid(f,u[0,i],u[1,i],h)
        else:
            print('no method was chosen')
            break
            
    return u
    
def euler(f,ti,yi,h):
    y = yi + h*f(ti, yi)
    
    return y

def heun(f,ti,yi,h):
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
    s1 = euler(f, ti, yi, h)
    y = yi + h/2*(f(ti, yi) + f(ti+h, s1))
    
    return y