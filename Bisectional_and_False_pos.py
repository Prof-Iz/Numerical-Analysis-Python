# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 03:42:04 2020

@author: Iz
"""

import pandas as pd #to use to store values
import numpy as np


def bisectional_method(f,xl,xu,iterations,valTrue=False):
    '''
    

    Parameters
    ----------
    f : Function of x - f(x)
        Pass in declared python function that takes in a single value.
    xl : Number
        Lower bound.
    xu : Number
        Upper Bound.
    iterations : integer
        How many iterations do you wish for the program to go through.
    valTrue : number, optional
        If Epsilon_T required providde value. The default is False.

    Returns
    -------
    Pandas Dataframe
        If successful returns a Pandas Data frame with all the important values.

    '''
    
    if (f(xu)*f(xl) > 0 ):
        return "Error"
    
    first = True  
    epsilon_a = [0]
    epsilon_t = []
    past = 0
    
    for i in range(iterations):
        
        xr = (xl + xu) / 2
        
        results = {"iteration" : i+1, 
                   "xl" : xl,
                   "xu" : xu,
                   "xr" : xr,
                   "F(xl)": f(xl),
                   "F(xu)": f(xu),
                   "F(xr)": f(xr),}
        if (first):
            resultsTable = pd.DataFrame(results,index=['iteration'])
        else:
            resultsTable = resultsTable.append(results,ignore_index=True)
        
        if (f(xr)*f(xl) < 0):
            xu = xr
        else:
            xl = xr
            
        if (not first):
            epsilon_a.append(abs((xr-past)/xr)*100)
            past = xr
        else:
            past = xr
            
        
        if (valTrue):
            epsilon_t.append(abs((valTrue-xr)/valTrue)*100)
        
        first = False
        
    
    resultsTable['epsilon_a'] = pd.Series(epsilon_a)
    if (valTrue):
        resultsTable['epsilon_t'] = pd.Series(epsilon_t)
        
    return resultsTable

def false_position(f,xl,xu,iterations,valTrue=False):
    '''
    

    Parameters
    ----------
    f : Function of x - f(x)
        Pass in declared python function that takes in a single value.
    xl : Number
        Lower bound.
    xu : Number
        Upper Bound.
    iterations : integer
        How many iterations do you wish for the program to go through.
    valTrue : number, optional
        If Epsilon_T required providde value. The default is False.

    Returns
    -------
    Pandas Dataframe
        If successful returns a Pandas Data frame with all the important values.

    '''
    
    if (f(xu)*f(xl) > 0 ):
        return "Error"
    
    first = True  
    epsilon_a = [0]
    epsilon_t = []
    past = 0
    
    for i in range(iterations):
        
        xr = xu - (f(xu)*(xl-xu))/(f(xl)-f(xu))
        
        results = {"iteration" : i+1, 
                   "xl" : xl,
                   "xu" : xu,
                   "xr" : xr,
                   "F(xl)": f(xl),
                   "F(xu)": f(xu),
                   "F(xr)": f(xr),}
        if (first):
            resultsTable = pd.DataFrame(results,index=['iteration'])
        else:
            resultsTable = resultsTable.append(results,ignore_index=True)
        
        if (f(xr)*f(xl) < 0):
            xu = xr
        else:
            xl = xr
            
        if (not first):
            epsilon_a.append(abs((xr-past)/xr)*100)
            past = xr
        else:
            past = xr
            
        
        if (valTrue):
            epsilon_t.append(abs((valTrue-xr)/valTrue)*100)
        
        first = False
        
    
    resultsTable['epsilon_a'] = pd.Series(epsilon_a)
    if (valTrue):
        resultsTable['epsilon_t'] = pd.Series(epsilon_t)
        
    return resultsTable
    
