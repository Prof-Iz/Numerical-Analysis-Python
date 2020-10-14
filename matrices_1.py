# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 02:49:18 2020

@author: Iz
"""

import numpy as np
from sympy.core.evaluate import evaluate
from sympy.matrices import Matrix,zeros,diag

def determinant(matrix):
    
    print("For matrix:\n")
    display(matrix)
    
    A = Matrix([[matrix[1,1],matrix[1,2]],[matrix[2,1],matrix[2,2]]])
    with evaluate(False):
        Aval = A[0,0]*A[1,1] - A[0,1]*A[1,0]
    print("The First Minor")
    display(A)
    display(Aval)
    display(Aval.doit())
    print("\n")
    
    B = Matrix([[matrix[1,0],matrix[1,2]],[matrix[2,0],matrix[2,2]]])
    with evaluate(False):
        Bval = B[0,0]*B[1,1] - B[0,1]*B[1,0]
    print("The Second Minor")
    display(B)
    display(Bval)
    display(Bval.doit())
    print("\n") 
    
    C = Matrix([[matrix[1,0],matrix[1,1]],[matrix[2,0],matrix[2,1]]])
    with evaluate(False):
        Cval = C[0,0]*C[1,1] - C[0,1]*C[1,0]
    print("The Third Minor")
    display(C)
    display(Cval)
    display(Cval.doit())
    print("\nDeterminant is ")
    det = matrix[0,0]*Aval.doit() - matrix[0,1]*Bval.doit() + matrix[0,2]*Cval.doit()
    display(det)
    return(det)
    
def cramer(m):
    
    det = determinant(m)
    
    x1 = Matrix([[m[0,3],m[0,1],m[0,2]],
                 [m[1,3],m[1,1],m[1,2]],
                 [m[2,3],m[2,1],m[2,2]]])
    
    print("x1 is :")
    
    
    display(x1)
    print("divided by det")
    display(det)
    print("is equal to")  
    display((x1.det()/det).doit())
    
    x2 = Matrix([[m[0,0],m[0,3],m[0,2]],
                 [m[1,0],m[1,3],m[1,2]],
                 [m[2,0],m[2,3],m[2,2]]])
    
    print("x2 is :")
    
    
    display(x2)
    print("divided by det")
    display(det)
    print("is equal to")  
    display((x2.det()/det).doit())
    
    x3 = Matrix([[m[0,0],m[0,1],m[0,3]],
                 [m[1,0],m[1,1],m[1,3]],
                 [m[2,0],m[2,1],m[2,3]]])
    
    print("x3 is :")
    
    
    display(x3)
    print("divided by det")
    display(det)
    print("is equal to")  
    display((x3.det()/det).doit())
    
def forward_elimination(m):
    print(f"R2 - {m[1,0]}/{m[0,0]} * R1 --> R2")
    try:
        m[1,:] = m[1,:] - (m[1,0]/m[0,0]) * m[0,:]
    except:
        null
    display(m)
    
    print(f"R3 - {m[2,0]}/{m[0,0]} * R1 --> R3")
    try:    
        m[2,:] = m[2,:] - (m[2,0]/m[0,0]) * m[0,:]
    except:
        null
    display(m)
    
    print(f"R3 - ({m[2,1]})/({m[1,1]}) * R2 --> R3")
    try:
        m[2,:] = m[2,:] - (m[2,1]/m[1,1]) * m[1,:]
    except:
        null
    display(m)
    
    try:
        x3 = m[2,3]/m[2,2]
    except:
        x3 = 0
        
    print(f"x3 = ({m[2,3]})/({m[2,2]}) = {x3}")
    
    try:
        x2 = ( m[1,3] - x3 * m[1,2]) / (m[1,1])
    except:
        x2 = 0
    
    print(f"x2 =  ({m[1,3]}) - ({x3}) * ({m[1,2]}) / ({m[1,1]})  = {x2}")
    
    try:
        x1 = (m[0,3] - x3*m[0,2] - x2*m[0,1]) / m[0,0]
    except:
        x1 = 0
    
    print(f"x1 = ({m[0,3]}) - ({x3}*{m[0,2]}) - ({x2}*{m[0,1]}) / ({m[0,0]}) = {x1}")
    
def LU(m):
    B = m[:,3]
    U = zeros(3,3)
    L = diag([1,1,1],unpack=True)
    U[0,:] = m[0,0:3]
    # Solving L and U
    print("For input Matrix:")
    display(m)
    
    l21 = m[1,0] / m[0,0]
    
    print(f"l21 = ({m[1,0]}) / ({m[0,0]}) = {l21}")
    
    l31 = m[2,0] / m[0,0]
    
    print(f"l31 = ({m[2,0]}) / ({m[0,0]}) = {l31}")
    
    U[1,:] = m[1,0:3] - l21 * m[0,0:3]
    U[2,:] = m[2,0:3] - l31 * m[0,0:3]
    
    
    print(f"R2 - {l21}*R1 --> R2")
    print(f"R3 - {l31}*R1 --> R3")
    print("U = ")
    display(U)
    
    l32 = U[2,1] / U[1,1]
    
    print(f"l32 = ({U[2,1]}) / ({U[1,1]})")
    print(f"R3 - {l32}*R2 --> R3")
    
    U[2,:] = U[2,:] - l32 * U[1,:]
    print("U = ")
    display(U)
    
    L[1,0] = l21
    L[2,0] = l31
    L[2,1] = l32
    
    print("L = ")
    display(L)
    
    # Solving intermediate D
    d1 = m[0,3]
    d2 = m[1,3] - L[1,0]*d1
    d3 = m[2,3] - L[2,0]*d1 - L[2,1]*d2
    
    D = Matrix([[d1],[d2],[d3]])
    print(f"""
    d1 = {m[0,3]}
    d2 = ({m[1,3]}) - ({L[1,0]})*({d1}) = {d2}
    d3 = ({m[2,3]}) - ({L[2,0]})*({d1}) - ({L[2,1]})*({d2}) = {d3}
    
    D = """)
    
    display(D)
    
    # Get X by back Substitution
    
    x3 = d3 / U[2,2]
    x2 = (d2 - x3* U[1,2]) / U[1,1]
    x1 = (d1 - x3*U[0,2] - x2*U[0,1]) / U[0,0]
    
    X = Matrix([[x1],[x2],[x3]])
    print(f"""
    x3 = ({d3}) / ({U[2,2]}) = {x3}
    x2 = (({d2}) - ({x3}* {U[1,2]}) / {U[1,1]} = {x2}
    x1 = (({d1}) - ({x3}*{U[0,2]}) - ({x2}*{U[0,1]})) / {U[0,0]} = {x1}
    X = 
    
    """)
    display(X)
    
def LU_inverse(m):
    B = m[:,3]
    U = zeros(3,3)
    L = diag([1,1,1],unpack=True)
    DA = L[:,0]
    DB = L[:,1] 
    DC = L[:,2]
    D_list = [DA,DB,DC]
    
    U[0,:] = m[0,0:3]
    # Solving L and U
    print("For input Matrix:")
    display(m)
    
    l21 = m[1,0] / m[0,0]
    
    print(f"l21 = ({m[1,0]}) / ({m[0,0]}) = {l21}")
    
    l31 = m[2,0] / m[0,0]
    
    print(f"l31 = ({m[2,0]}) / ({m[0,0]}) = {l31}")
    
    U[1,:] = m[1,0:3] - l21 * m[0,0:3]
    U[2,:] = m[2,0:3] - l31 * m[0,0:3]
    
    
    print(f"R2 - {l21}*R1 --> R2")
    print(f"R3 - {l31}*R1 --> R3")
    print("U = ")
    display(U)
    
    l32 = U[2,1] / U[1,1]
    
    print(f"l32 = ({U[2,1]}) / ({U[1,1]})")
    print(f"R3 - {l32}*R2 --> R3")
    
    U[2,:] = U[2,:] - l32 * U[1,:]
    print("U = ")
    display(U)
    
    L[1,0] = l21
    L[2,0] = l31
    L[2,1] = l32
    
    print("L = ")
    display(L)
    Xn = []
    # Solving intermediate D
    for D_i in D_list:
        print("For :  ")
        display(D_i)
        
        d1 = D_i[0]
        d2 = D_i[1] - L[1,0]*d1
        d3 = D_i[2] - L[2,0]*d1 - L[2,1]*d2
        
        D = Matrix([[d1],[d2],[d3]])
        print(f"""
        d1 = {D_i[0]}
        d2 = ({D_i[1]}) - ({L[1,0]})*({d1}) = {d2}
        d3 = ({D_i[2]}) - ({L[2,0]})*({d1}) - ({L[2,1]})*({d2}) = {d3}
        
        D = """)
        
        display(D)
        
        # Get X by back Substitution
        
        x3 = d3 / U[2,2]
        x2 = (d2 - x3* U[1,2]) / U[1,1]
        x1 = (d1 - x3*U[0,2] - x2*U[0,1]) / U[0,0]
        
        X = Matrix([[x1],[x2],[x3]])
        Xn.append(X)
        print(f"""
        x3 = ({d3}) / ({U[2,2]}) = {x3}
        x2 = (({d2}) - ({x3}* {U[1,2]}) / {U[1,1]} = {x2}
        x1 = (({d1}) - ({x3}*{U[0,2]}) - ({x2}*{U[0,1]})) / {U[0,0]} = {x1}
        X = 
        
        """)
        display(X)
    
    print("Inverse Matrix is ")
    Xf = zeros(3,3)
    i = 0
    for col in Xn:
        Xf[:,i] = col
        i+=1
    display(Xf)
        
    