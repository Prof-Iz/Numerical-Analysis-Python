# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 09:46:06 2020

@author: Iz
"""

import pandas as pd

def secant_method(f,xi,xi_minus1,iterations,valTrue=False):
    
    '''
    Function to find derivative using  Secant method

    Parameters
    ----------
    f : function of x defined in python def f(x):

    xi : number
        Initial point.
        
    d : number
        Step.
    iterations : number
        how many iterations.
    valTrue : number, optional
        Enter True value if Epsilon_t is 
        required. The default is False.

    Returns
    -------
    resultsTable : Pandas DataFrame
        Table containing all the results.

    '''
    
    first = True  
    epsilon_a = [0]
    epsilon_t = []
    past = 0
    
    for i in range(iterations):
        
        function_xi = f(xi)
        function_xi_minus1 = f(xi_minus1)
        
        results = {"iteration" : i+1, 
                   "xi" : xi,
                   "f(xi)" : function_xi,
                   r"x_{i-1}" : xi_minus1,
                   r"f(x_{i-1})": function_xi_minus1,}
        

        past = xi
        xNext = xi - (function_xi*(xi_minus1-xi))/(function_xi_minus1-function_xi)
        xi , xi_minus1 = xNext , xi
        
        
        
        
        if (first):
            resultsTable = pd.DataFrame(results,index=['iteration'])
        else:
            resultsTable = resultsTable.append(results,ignore_index=True)
            epsilon_a.append(abs((xi-past)/xi)*100)

        if (valTrue):
            epsilon_t.append(abs((valTrue-xi)/valTrue)*100)
        
        first = False
        
    resultsTable['epsilon_a'] = pd.Series(epsilon_a)
    if (valTrue):
        resultsTable['epsilon_t'] = pd.Series(epsilon_t)
        
    return resultsTable
        
    


def mod_secant_method(f,xi,d,iterations,valTrue=False):
    '''
    Function to find derivative using modified Secant method

    Parameters
    ----------
    f : function of x defined in python def f(x):

    xi : number
        Initial point.
        
    d : number
        Step.
    iterations : number
        how many iterations.
    valTrue : number, optional
        Enter True value if Epsilon_t is 
        required. The default is False.

    Returns
    -------
    resultsTable : Pandas DataFrame
        Table containing all the results.

    '''
    
    first = True  
    epsilon_a = [0]
    epsilon_t = []
    past = 0
    
    for i in range(iterations):
        
        
        
        results = {"iteration" : i+1, 
                   "xi" : xi,
                   "f(xi)" : f(xi),
                   "xi+d" : xi+d,
                   "f(xi + d)": f(xi + d),}
        

        past = xi
        xi = xi-((d*xi*f(xi))/(f(xi+(d*xi))-f(xi)))
        
        
        if (first):
            resultsTable = pd.DataFrame(results,index=['iteration'])
        else:
            resultsTable = resultsTable.append(results,ignore_index=True)
            epsilon_a.append(abs((xi-past)/xi)*100)

        if (valTrue):
            epsilon_t.append(abs((valTrue-xi)/valTrue)*100)
        
        first = False
        
    resultsTable['epsilon_a'] = pd.Series(epsilon_a)
    if (valTrue):
        resultsTable['epsilon_t'] = pd.Series(epsilon_t)
        
    return resultsTable
        