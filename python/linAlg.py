import numpy as np

def vectorNorm(x,case):
    if(case == 2):
        norm = np.sqrt(np.sum([abs(i)**2 for i in x]))
        return norm
    else:
        norm = np.array([abs(i) for i in x])
        if(case == 0):
            return norm.max()
        elif(case == 1):
            return norm.sum()
        

def jacobi(A,b,x0,tol):
    U = np.triu(A,1)
    L = np.tril(A,-1)
    D = np.diag(np.diag(A))
    
    k = 0
    x = np.linalg.solve(-D,np.matmul(L+U,x0))
    x += np.linalg.solve(D,b)
        
    r = []
    r.append(np.linalg.norm(b-np.matmul(A,x)))
    cond = True
    
    while(r[k] >= tol):
        x = np.linalg.solve(-D,np.matmul(L+U,x))
        x += np.linalg.solve(D,b)
        r.append(np.linalg.norm(b-np.matmul(A,x)))
        k += 1
        
    return x,r,k

def gauss_seidel(A,b,x,tol):
    L = np.tril(A)
    U = A - L
    
    x = np.dot(np.linalg.inv(L), b - np.dot(U, x))
    
    k = 0
    r = []
    r.append(np.linalg.norm(b-np.matmul(A,x)))
    while(r[k] >= tol):
        x = np.dot(np.linalg.inv(L), b - np.dot(U, x))
        r.append(np.linalg.norm(b-np.matmul(A,x)))
        k += 1
    return x,r,k

def lstSqVals(x,y,basis):
    A = np.array(list(map(basis,x)))
    
    rows,cols = A.shape
    
    vals = np.linalg.lstsq(A,y,rcond=None)
    c = vals[0]
    rms = vals[1]/np.sqrt(rows)
    
    return c,rms[0]

def runLstSqVals(x,p,t,f,base):
    f = fexact(t)

    for i in range(len(x)):
    g = 0
    for j in range(5):
        y = fexact(x[i])
        A = base(x[i],j)
        
        vals = np.linalg.lstsq(A,y,rcond=None)
        c = vals[0]
        rms = vals[1]/np.linalg.norm(y)
        
        p = np.polyval(c[::-1],t)        
        
        res = np.linalg.norm(f-p)/np.linalg.norm(f)
        msg = 'M = {0}, P = {1}: (1) RES {2:1.3e}, (2) RES {3:1.3e}'.format(len(x[i]),j+1,rms[0],res)
        print(msg)