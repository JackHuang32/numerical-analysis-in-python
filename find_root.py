import math
from math import tan,pi,log,exp
from sympy import*
def function1(x):
    return (x/(1+x**2))-(500/841)*(1-21*x/125)
def function2(x):
    return x**4-x**3-3*(x**2)+5*x-2
def bis(f,init_a,init_b):
    a = init_a;b = init_b
    p = (a+b)/2
    arr = []
    for i in range(10):
        if f(p)*f(a) < 0:
            b = p
        else:
            a = p
        p = (a+b)/2
        arr.append(p)
    return arr
def d(f,v):
    f = f.diff(x)
    der = lambdify(x,f)
    return der(v)
def newton(f,init,times):
    arr=[]
    cur = init
    func = lambdify(x,f)
    for i in range(times):
        cur = cur-func(cur)/d(f,cur)
        arr.append(cur)
    return arr
def secant(f,init_a,init_b):
    p1 = init_a;p2 = init_b
    arr = []
    for i in range(5):
        p = p2 - f(p2)*(p2-p1)/(f(p2)-f(p1))
        p1 = p2
        p2 = p
        arr.append(p2)
    return arr
x = symbols('x')
f = (x**3)*(1-0.5**2)+(0.4*4.2*(0.5**2)-150*(4.2**2)/(40**2))*(x**2)\
    +(150**2)*(4.2**4)*x/(3*(40**4))-((150*(4.2**2))/3*(40**2))**3




