# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 16:20:47 2020

@author: Iz
"""
import pandas as pd
import sympy as sp

x = sp.symbols("x")

def fixed_point(f,xo,iterations,valTrue = False):
    '''
    

    Parameters
    ----------
    f : Function in terms of X
        
    xo : number
        Beginning point for iteration.
    iterations : number
        how many iterations required.
    valTrue : number, optional
        Replace with true value if Epsilon
        t is required. The default is False.

    Returns
    -------
    resultsTable : Pandas dataFrame
        table of results.

    '''
    
    first = True  
    epsilon_a = [0]
    epsilon_t = []
    past =0
    
    for i in range(iterations):
        
        results = {"iteration":i+1,
                   "xi":xo,
                   "f(xi) = xi+1":f(xo)}
        
        past = xo
        xo = f(xo)
        
        if (first):
            resultsTable = pd.DataFrame(results,index=['iteration'])
        else:
            resultsTable = resultsTable.append(results,ignore_index=True)
            epsilon_a.append(abs((xo-past)/xo)*100)

        if (valTrue):
            epsilon_t.append(abs((valTrue-xo)/valTrue)*100)
        
        first = False
        
    resultsTable['epsilon_a'] = pd.Series(epsilon_a)
    if (valTrue):
        resultsTable['epsilon_t'] = pd.Series(epsilon_t)
        
    return resultsTable

def newton_raphson(f,xo,iterations,valTrue=False,precision=5):
    '''
    

    Parameters
    ----------
    f : f = equation in terms of sympy symbols
        
    xo : number
        Beginning point for iteration.
        
    iterations : number
        how many iterations required.
        
    valTrue : number, optional
        Replace with true value if Epsilon
        t is required. The default is False.
        
    precision : int, optional
        precision of decimal places. The default is 5.

    Returns
    -------
    resultsTable : Pandas Dataframe
        Datafram containing results.

    '''
    
    first = True  
    epsilon_a = [0]
    epsilon_t = []
    past =0
    
    for i in range(iterations):
        
        func = f.subs(x,xo).evalf(5)
        func_prime = sp.diff(f,x).subs(x,xo).evalf(precision)
        
        results = {"iteration":i+1,
                   "xi":xo,
                   "f(xi)":func,
                   "f'(xi)":func_prime}
        
        past = xo
        xo = xo - (func/func_prime)
        
        
        if (first):
            resultsTable = pd.DataFrame(results,index=['iteration'])
        else:
            resultsTable = resultsTable.append(results,ignore_index=True)
            epsilon_a.append(abs((xo-past)/xo)*100)

        if (valTrue):
            epsilon_t.append(abs((valTrue-xo)/valTrue)*100)
        
        first = False
        
    resultsTable['epsilon_a'] = pd.Series(epsilon_a)
    if (valTrue):
        resultsTable['epsilon_t'] = pd.Series(epsilon_t)
        
    return resultsTable