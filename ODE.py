import numpy as np
from math import factorial
from sympy import*
def Euler(f,a,b,n,init):
    t_ = np.linspace(a,b,n+1)
    y=[init]
    for i in range(1,n+1):
        y.append(y[i-1]+((b-a)/n)*lambdify((t,x),f)(t_[i-1],y[i-1]))
    return y
def d(f,n,t_,x_):
    for i in range(n):
        f = diff(f,t)+diff(f,x)*f
    return lambdify((t,x),f)(t_,x_)
def Taylor(f,a,b,n,init,order):
    h = (b-a)/n
    t_ = np.linspace(a,b,n+1)
    y = [init]
    for i in range(1,n+1):
        tmp = 0
        for j in range(1,order+1):
            tmp+=(h**j*d(f,j-1,t_[i-1],y[i-1]))/factorial(j)
        y.append(y[i-1]+tmp)
    return y
def RK(M,f_,n,a,b,init):
    f = lambdify((t,x),f_)
    w=[init]
    h = (b-a)/n
    domain = np.linspace(a,b,n+1)
    for i in range(n):
        k=[]
        for j in range(M.shape[0]-1):
            ak=0
            for l in range(j):
                ak+=M[j][l+1]*k[l]
            k.append(h*f(domain[i]+M[j][0]*h,w[i]+ak))
        w.append(w[i]+sum([M[-1][a+1]*k[a] for a in range(M.shape[0]-1)])) 
    return w    
x = symbols('x')
t = symbols('t')
f = 1+x/t
M = np.array([[0,0,0],[1,1,0],[0,1/2,1/2]],dtype=float)
print(RK(M,f,10,1,6,1))
