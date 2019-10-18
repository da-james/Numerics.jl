import numpy as np

def euler(f,ti,tf,yi,n=0,h=0):
    if(h == 0 and n != 0):
        h = (tf-ti)/n
    elif(h != 0 and n == 0):
        n = int((tf-ti)/h)
        
    y = np.zeros(n+1)
    t = np.zeros(n+1)
    
    y[0] = yi
    t[0] = ti
    
    for i in range(n):
        t[i+1] = t[i] + h
        y[i+1] = y[i] + h*f(t[i],y[i])
        
    return t,y

def heun(f,ti,tf,yi,n=0,h=0):
    if(h == 0 and n != 0):
        h = (tf-ti)/n
    elif(h != 0 and n == 0):
        n = int((tf-ti)/h)
    
    y = np.zeros(n+1)
    t = np.zeros(n+1)
    
    y[0] = yi
    t[0] = ti
    
    for i in range(n):
        s1 = f(t[i],y[i])
        s2 = f(t[i] + h,y[i] + h*s1)
        
        t[i+1] = t[i] + h
        y[i+1] = y[i] + h/2*(s1 + s2)
        
    return t,y