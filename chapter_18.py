import numpy as np

def linear_interpolation(f,x0,x1,xn,valTrue=False):
    f1 = f(x0) + (f(x1)-f(x0))/(x1-x0) * (xn-x0)
    print(f'''
    f1 = f(x0) + (f(x1)-f(x0))/(x1-x0) * (x-x0)
    = {f(x0)} + ({f(x1)}-{f(x0)})/({x1-x0}) * ({xn-x0}) = {f1}
    ''')

    if (valTrue):
        epsilon_t = abs((valTrue - f1) / valTrue * 100)
        print(f'True Value Error is {epsilon_t}%')


def quadratic_interpolation(x0,x1,x2,f0,f1,f2,x,valTrue=False):
    b0 = f0
    b1 = (f1 -f0)/(x1-x0)
    b2 = ((f2 - f1)/(x2-x1) - b1)/(x2-x0)
    def g(x):
        return b0 + b1*(x-x0)+b2*(x-x0)*(x-x1)

    answer = g(x)
    print(f"""
    b0 = f0 = {b0}
    b1 = (f1 -f0)/(x1-x0) = ({f1} -{f0})/({x1-x0})
    = {b1}

    b2 = ((f2 - f1)/(x2-x1) - b1)/(x2-x0)
    = (({f2 - f1})/({x2-x1}) - {b1})/({x2-x0})
    = {b2}
    def g(x):
        return b0 + b1*(x-x0)+b2*(x-x0)(x-x1)
    
    {answer}
    
    """)

    if (valTrue):
        epsilon_t = abs((valTrue - answer ) / valTrue * 100)
        print(f'True Value Error is {epsilon_t}%')