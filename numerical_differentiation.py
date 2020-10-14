# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 03:00:25 2020

@author: Izdhan
"""

def forward_finite(f,xi,h,trueVal=False):
    '''
    The Forward finite method of differentiating an equation
    Numerically.

    Parameters
    ----------
    f : function
        Declare and pass your f(x) here.
    xi : number
        Point about which you wish to check gradient.
    h : number
        Step Size between the methods.
    trueVal : number , optional
        True value of calculation if present. The default is False.

    Returns
    -------
    f_prime : number
        numerical value of f'(x) about provided point.

    '''
    xn = xi + h
    f_prime = (f(xn)-f(xi))/(xn-xi)
    
    if (trueVal):
        epsilon_t = abs((trueVal - f_prime) / trueVal * 100)
        print(f'True Value Error is {epsilon_t}%')
        
    return f_prime

def backward_finite(f,xi,h,trueVal=False):
    '''
    The Backward finite method of differentiating an equation
    Numerically.

    Parameters
    ----------
    f : function
        Declare and pass your f(x) here.
    xi : number
        Point about which you wish to check gradient.
    h : number
        Step Size between the methods.
    trueVal : number , optional
        True value of calculation if present. The default is False.

    Returns
    -------
    f_prime : number
        numerical value of f'(x) about provided point.

    '''
    xb = xi - h
    f_prime = (f(xi)-f(xb))/(xi-xb)
    
    if (trueVal):
        epsilon_t = abs((trueVal - f_prime) / trueVal * 100)
        print(f'True Value Error is {epsilon_t}%')
        
        
    return f_prime

def center_finite(f,xi,h,trueVal=False):
    '''
    The Center finite method of differentiating an equation
    Numerically.

    Parameters
    ----------
    f : function
        Declare and pass your f(x) here.
    xi : number
        Point about which you wish to check gradient.
    h : number
        Step Size between the methods.
    trueVal : number , optional
        True value of calculation if present. The default is False.

    Returns
    -------
    f_prime : number
        numerical value of f'(x) about provided point.

    '''
    xn = xi + h
    xb = xi - h
    f_prime = (f(xn)-f(xb))/(xn-xb)
    
    if (trueVal):
        epsilon_t = abs((trueVal - f_prime) / trueVal * 100)
        print(f'True Value Error is {epsilon_t}%')
    
    return f_prime
    
def numerical_diff(f,xi,h,trueVal=False):
    '''
    
    Combines all three numerical methods of differentiation
    Forward finite - Center finite - Backward Finite
    and performs all three of these methods.
    
    If you wish to get a value returned use the seperate
    function calls for each method
    
    forward_finite()
    center_finite()
    backward_finite()

    Parameters
    ----------
    f : function
        Declare and pass your f(x) here.
    xi : number
        Point about which you wish to check gradient.
    h : number
        Step Size between the methods.
    trueVal : number , optional
        True value of calculation if present. The default is False.

    Returns
    -------
    None.

    '''
    
    forward = forward_finite(f,xi,h,trueVal)
    print(f"""Value of f'({xi}) = {forward}
          Result above by forward finite method\n""")
          
    backward = backward_finite(f,xi,h,trueVal)
    print(f"""Value of f'({xi}) = {backward}
          Result above by backward finite method\n""")
          
    center = center_finite(f,xi,h,trueVal)
    print(f"""Value of f'({xi}) = {center}
          Result above by center finite method\n""")