# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 21:12:25 2020

@author: Izdhan
"""


import sympy as sp
import numpy as np

x, y = sp.symbols("x y")

# insert function
g = sp.Eq(x**2*y-1,0)
display(g)
f = sp.Derivative(g,x)
print(f)

def taylors_method(f,xn,xo,yo,derivative=False):
    h = xn - xo
    
    if (not derivative):
        fxn = yo + sp.diff(f, x).subs(x,xo)*h + sp.diff(f,x,2).subs(x,xo)*(h**2/np.math.factorial(2)) + sp.diff(f,x,3).subs(x,xo)*(h**3/np.math.factorial(3)) + sp.diff(f,x,4).subs(x,xo)*(h**4/np.math.factorial(4))
        print(fxn)
    