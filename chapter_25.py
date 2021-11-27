# -*- coding: utf-8 -*-
"""
Created on Sun Nov 29 16:36:43 2020

@author: User
"""

def runge_kutta_4(f,x0,y0,h,n):
    for i in range(0,n):
        k1 = h*f(x0,y0)
        k2 = h*f(x0+h/2,y0+k1/2)
        k3 = h*f(x0+h/2,y0+k2/2)
        k4 = h*f(x0+h,y0+k3)
        
        dy = 1/6 * (k1 + 2*k2 +2*k3 + k4)
        
        yn = y0 + dy
        
        print(f"""
        xi = {x0}, yi = {y0}      
        
        k1 = {h}*f(x0,y0) = {h}*{f(x0,y0)} = {k1}
        k2 = h*f(x0+h/2,y0+k1/2) = {h}*{f(x0+h/2,y0+k1/2)} = {k2}
        k3 = h*f(x0+h/2,y0+k2/2) = {h}*{f(x0+h/2,y0+k2/2)} = {k3}
        k4 = h*f(x0+h,y0+k3) = {h}*{f(x0+h,y0+k3)} = {k4}
        
        dy = 1/6 * (k1 + 2*k2 +2*k3 + k4) = 1/6 * ({k1} + {2*k2} +{2*k3} + {k4})
        dy = {dy}
        
        yn = y0 + dy = {y0} + {dy}
        
        """)
        
        print(f"x = {x0+h}, y = {yn}")
        
        x0 = x0+h
        y0 = yn
        

def euler_ode(fxy,x0,y0,xn,n=1):
    '''
    

    Parameters
    ----------
    fxy : function (x,y)
        Python function that takes x,y and returns
        a scalar value.
    x0 : number
        initial X
    y0 : number
        initial y.
    xn : number
        x to predict, this will set the h for subsequent
        iterations.
    n : number
        Number of iterations.Default is 1 iteration


    Returns
    -------
    yn : number
        The predicted y-value after n increments.

    '''
    
    for i in range(0,n):
    
        FXY = fxy(x0,y0)
        h = xn - x0
        
        yn = y0 + FXY*h
        
        print(f"""
              
        FXY = fxy(x0,y0) = {FXY}
        h = xn - x0 = {h}
        
        yn = y0 + FXY*h = {y0} + {FXY}*{h}
        yn = {yn}    
        
        
        """)    
   
        x0 , xn, y0 = xn, xn+h, yn    
    return yn
        
    

def heun(fxy,x0,y0,h,n=1):
    for i in range(0,n):
        
        FXY = fxy(x0,y0)
        y01 = y0 + FXY*h
        
        xn = x0 +h
        
        FXY1 = fxy(xn,y01)
        
        yn = y0 + h*((FXY+FXY1)/2)
        
        print(f"""
         
        FXY = f(x0,y0) = {FXY}
        y{i+1} = y{i} + FXY*h = {y01}
        
        x{i+1} = x0 +h = {xn}
        
        FXY{i+1} = f(xn,y{i+1}) = {FXY1}
        
        yn = y{i} + h*((FXY+FXY1)/2) = {y0} + h*(({FXY}+{FXY1})/2) = {yn}
    
        
        """)
        
        x0 , y0 = xn, yn
        

def runge_kutta_5(f,x0,y0,h,n=1):
    for i in range(0,n):
        k1 = f(x0,y0)
        k2 = f(x0+h/4,y0 + 1/4*k1*h)
        k3 = f(x0+h/4,y0+1/8*k1*h+1/8*k2*h)
        k4 = f(x0+h/2,y0-1/2*k2*h+k3*h)
        k5 = f(x0+3*h/4,y0+3/16*k1*h+9/16*k4*h)
        k6 = f(x0+h,y0-3/7*k1*h+2/7*k2*h+12/7*k3*h-12/7*k4*h+8/7*k5*h)
        
        
        dy = 1/90 * (7*k1 +32*k3 + 12*k4+32*k5 + 7*k6) * h
        
        yn = y0 + dy
        
        print(f"""
        xi = {x0}, yi = {y0}      
        
        k1 = f(x0,y0) ={k1}
        k2 = f({x0+h/4},{y0 + 1/4*k1*h}) = {k2}
        k3 = f(x0+h/4,y0+1/8*k1*h+1/8*k2*h) = {k3}
        k4 = f({x0+h/2},{y0-1/2*k2*h+k3*h}) = {k4}
        k5 = f({x0+3*h/4},{y0+3/16*k1*h+9/16*k4*h}) = {k5}
        k6 = f({x0+h},{y0-3/7*k1*h+2/7*k2*h+12/7*k3*h-12/7*k4*h+8/7*k5*h}) = {k6}
        
        
        dy = 1/90 * ({7*k1 +32*k3 + 12*k4+32*k5 + 7*k6}) * {h}
        dy = {dy}
        yn = y{i} + dy = {yn}
      
        """)
        
        print(f"x = {x0+h}, y = {yn}")
        
        x0 = x0+h
        y0 = yn
        
        
def runge_kutta_3(f,x0,y0,h,n=1,val_True=False):
    for i in range(0,n):
        k1 = f(x0,y0)
        k2 = f(x0+h/2,y0 + 1/2*k1*h)
        k3 = f(x0+h,y0-k1*h+2*k2*h)
    
        
        
        dy = 1/6 * (k1 + 4*k2 + k3) * h
        
        yn = y0 + dy
        
        print(f"""
        xi = {x0}, yi = {y0}      
        
        k1 = f(x0,y0)
        k2 = f({x0+h/2},{y0 + 1/2*k1*h}) = {k2}
        k3 = f({x0+h},{y0-k1*h+2*k2*h}) = {k3}
    
        
        
        dy = 1/6 * (k1 + 4*k2 + k3) * h = {dy}
        
        yn = y0 + dy = {yn}
      
        """)
        
        print(f"x = {x0+h}, y = {yn}")
        
        x0 = x0+h
        y0 = yn
        
        
def runge_kutta_2_midpoint(f,x0,y0,h,n=1):
    
    for i in range(0,n):
        k1 = f(x0,y0)
        k2 = f(x0+h/2,y0 + 1/2*k1*h)
    
        
        
        dy = k2*h
        
        yn = y0 + dy
        
        print(f"""
        xi = {x0}, yi = {y0}      
        
        k1 = f(x0,y0) = {k1}
        k2 = f({x0+h/2},{y0 + 1/2*k1*h}) = {k2}
    
        dy = k2*h = {dy}
        
        yn = y0 + dy = {yn}
      
        """)
        
        print(f"x = {x0+h}, y = {yn}")
        
        x0 = x0+h
        y0 = yn
        
def runge_kutta_2_ralston(f,x0,y0,h,n=1):
    for i in range(0,n):
        k1 = f(x0,y0)
        k2 = f(x0+3*h/4,y0 + 3/4*k1*h)
    
        
        
        dy = (1/3*k1+2/3*k2)*h
        
        yn = y0 + dy
        
        print(f"""
        xi = {x0}, yi = {y0}      
        
        k1 = f(x0,y0) = {k1}
        k2 = f({x0+3*h/4},{y0 + 3/4*k1*h}) = {k2}
        
        dy = ({1/3*k1+2/3*k2})*h = {dy}
        
        yn = y0 + dy = {yn}
      
        """)
        
        print(f"x = {x0+h}, y = {yn}")
        
        x0 = x0+h
        y0 = yn