from  scipy.interpolate import CubicSpline
import numpy as np
from math import cos,pi,sin
#if you want P_0,1,2,3(value) with x,y then call Neville(x,y,value,[0,1,2,3])
#index start with 0 
#the function return P(value)
def Neville(x,y,value,arr):
    M = np.zeros((len(x[0]),len(x[0])),dtype=float)
    M[0] = y[0]
    for i in range(1,M.shape[0]):
        for j in range(i,M.shape[1]):
            entry = ((value-x[0][j-i])*M[i-1][j]-(value-x[0][j])*M[i-1][j-1])/(x[0][j]-x[0][j-i])
            M[i][j] = entry
    print(M.T)
    print('P',arr,'is',M[len(arr)-1,arr[0]+len(arr)-1])
    return M[M.shape[0]-1][M.shape[1]-1]
#if you want f[x1,x3,x5,x7], then call print(divided_diff(x,y,[1,3,5,7])
#index start with 0
def divided_diff(x,y,arr):
    if len(arr)==1:
        return y[0][arr[0]]
    else:
        return (divided_diff(x,y,arr[1:])-divided_diff(x,y,arr[:-1]))/(x[0][arr[-1]]-x[0][arr[0]])
def solve_tri(tri,b_,sol):
    P=[tri[i][i] for i in range(tri.shape[0])]
    Q=[tri[i][i+1] for i in range(tri.shape[0]-1)]
    R=[0 for i in range(len(P))]
    for i in range(len(Q)):
        R[i+1] = tri[i+1][i]
    Y=[0 for i in range(len(Q))]
    Y[0]=Q[0]/P[0]
    for i in range(1,len(Y)):
        Y[i]=Q[i]/(P[i]-R[i]*Y[i-1])
    W = [0 for i in range(len(P))]
    W[0]=b_[0]/P[0]
    for i in range(1,len(P)):
        W[i]=(b_[i]-R[i]*W[i-1])/(P[i]-R[i]*Y[i-1])
    sol[-1] = W[-1]
    for i in range(len(sol)-2,-1,-1):
        sol[i] = W[i]-Y[i]*sol[i+1]
#return a matrix of cubic spline coefficient with M[i][0]=ai M[i][1]=bi M[i][2]=ci M[i][3]=di
def Not_a_knot(x,y):
    #solve tridiagonal
    h = []
    for i in range(len(x[0])-1):
        h.append(x[0][i+1]-x[0][i])
    b_ = []
    for i in range(0,len(x[0])-2):
        b_.append(((3/h[i+1])*(y[0][i+2]-y[0][i+1])-(3/h[i])*(y[0][i+1]-y[0][i]))) 
    tri = np.zeros((len(x[0])-2,len(x[0])-2),dtype=float)
    sol = [0 for i in range(tri.shape[0])]#solve c1~cn-1
    tri[0][0] = 3*h[0]+2*h[1]+h[0]**2/h[1]
    tri[0][1] = h[1]-h[0]**2/h[1]
    tri[tri.shape[0]-1][tri.shape[1]-2] = h[-2]-h[-1]**2/h[-2]
    tri[tri.shape[0]-1][tri.shape[1]-1] = 3*h[-1]+2*h[-2]+h[-1]**2/h[-2]
    for i in range(1,tri.shape[0]-1):
        tri[i][i-1] = h[i]
        tri[i][i] = 2*(h[i]+h[i+1])
        tri[i][i+1] = h[i+1]
    solve_tri(tri,b_,sol)
    c = [];b=[];a=[];d=[0 for i in range(len(y[0])-1)]
    c.append((1+h[0]/h[1])*sol[0]-(h[0]/h[1])*sol[1])
    for i in sol:
        c.append(i)
    c.append(-(h[-1]/h[-2])*c[-2]+(1+h[-1]/h[-2])*c[-1])
    for i in range(len(c)-1):
        b.append((y[0][i+1]-y[0][i])/h[i]-(2*c[i]+c[i+1])*h[i]/3)
    for i in range(len(y[0])-1):
        a.append(y[0][i])
    for i in range(1,len(d)-1):
        d[i] = (c[i+1]-c[i])/(3*h[i])
    d[0] = d[1];d[-1] = d[-2]
    M = np.zeros((4,len(b)),dtype=float)
    M[0] = y[0][:-1];M[1] = np.array([b]);M[2] = np.array([c[:-1]]);M[3] = np.array([d])
    return M.T
#return s(value) where s is the cubic spline polynomial
def cal_not_a_knot(x,y,value):
    M = Not_a_knot(x,y)
    xi=0
    for i in range(len(x[0])-1):
        if value>=x[0][i] and value<x[0][i+1]:
            xi = i
    return M[xi][0]+M[xi][1]*(value-x[0][xi])+M[xi][2]*(value-x[0][xi])**2+M[xi][3]*(value-x[0][xi])**3    
#similar to Not_a_knot
def Clamped(x,y,der_a,der_b):
    h = []
    for i in range(len(x[0])-1):
        h.append(x[0][i+1]-x[0][i])
    b_=[]
    b_.append((3/h[0])*(y[0][1]-y[0][0])-3*der_a)
    for i in range(0,len(h)-1):
        b_.append(((3/h[i+1])*(y[0][i+2]-y[0][i+1])-(3/h[i])*(y[0][i+1]-y[0][i])))
    b_.append(3*der_b-(3/h[-1])*(y[0][-1]-y[0][-2]))
    tri = np.zeros((len(x[0]),len(x[0])),dtype=float)
    sol = [0 for i in range(tri.shape[0])]#solve c1~cn-1
    tri[0][0] = 2*h[0]
    tri[0][1] = h[0]
    tri[tri.shape[0]-1][tri.shape[1]-2] = h[-1]
    tri[tri.shape[0]-1][tri.shape[1]-1] = 2*h[-1]
    for i in range(1,tri.shape[0]-1):
        tri[i][i-1] = h[i-1]
        tri[i][i] = 2*(h[i-1]+h[i])
        tri[i][i+1] = h[i]
    solve_tri(tri,b_,sol)
    c = [];b=[];a=[];d=[0 for i in range(len(y[0])-1)]
    for i in sol:
        c.append(i)
    for i in range(len(c)-1):
        b.append((y[0][i+1]-y[0][i])/h[i]-(2*c[i]+c[i+1])*h[i]/3)
    for i in range(len(y[0])-1):
        a.append(y[0][i])
    for i in range(len(d)):
        d[i] = (c[i+1]-c[i])/(3*h[i])
    M = np.zeros((4,len(b)),dtype=float)
    M[0] = y[0][:-1];M[1] = np.array([b]);M[2] = np.array([c[:-1]]);M[3] = np.array([d])
    return M.T
def cal_clamped(x,y,der_a,der_b,value):
    M = Clamped(x,y,der_a,der_b)
    xi=0
    for i in range(len(x[0])-1):
        if value>=x[0][i] and value<x[0][i+1]:
            xi = i
    return M[xi][0]+M[xi][1]*(value-x[0][xi])+M[xi][2]*(value-x[0][xi])**2+M[xi][3]*(value-x[0][xi])**3    
x = np.array([[-1,-0.5,0,0.5,1]],dtype=float)
y = np.array([[0,0.82436,1,0.90980,0.73576]],dtype=float)

print(cal_clamped(x,y,2.71828,-0.36788,0.3))


