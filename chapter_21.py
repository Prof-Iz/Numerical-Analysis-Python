import numpy as np
import pandas as pd

def trapezoidal_rule(f,a,b,val_true=False,suppress_working=False):
    """Single trapezoidal Rule implementation

    Args:
        f (function declared in python f(x))
        a ([type]): [description]
        b ([type]): [description]
        val_true (bool, optional): [description]. Defaults to False.
    """
    
    fa = round(f(a),5)
    print(f"f(a) = {fa}")

    fb = round(f(b),5)
    print(f"f(b) = {fb}")

    I = (b-a)*((fa+fb)/2)
    I = round(I,5)
    print(f"I = ({b}-{a})*({fa})+{fb})/2")
    print(f"I = {I}")

    if val_true:
        epsilon_t = ((val_true - I)/ val_true)*100
        epsilon_t = round(epsilon_t,5)
        print(f"epsilon_t = ({val_true} - {I})/ {val_true})*100")
        print(f"e_t is {epsilon_t}")
   
    
    
def multi_trapezoidal_rule(f,a,b,n,val_true=False):
    '''
    

    Parameters
    ----------
    f : python function, f(x)
    a : number, starting x
    b : number, ending x
        
    n : number of iterations
        
    val_true : number, true value of area
        The default is False.

    Returns
    -------
    None.

    '''
    
    i=0
    h = (b-a)/n
    print(f"h = ({b}-{a})/{n}")
    working = pd.DataFrame({"x":a,"f(x)":f(a)},index=[i])
    for i in range(1,n+1):
        x = a+h*i
        working = working.append({"x":x,"f(x)":round(f(x),5)},ignore_index=True)
    print(working)
    sum_f =np.sum(working['f(x)'][1:-1])
    print(f"sum of that part is {sum_f}")
    I = (b-a)*((working['f(x)'][0]+2*sum_f+working['f(x)'][n])/(2*n))
    print(f"I = ({b}-{a})*(({working['f(x)'][0]}+2*{sum_f}+{working['f(x)'][n]})/(2*{n}))")
    print(f"I = {round(I,5)}")    
    
    if val_true:
        epsilon_t = ((val_true - I)/ val_true)*100
        epsilon_t = round(epsilon_t,5)
        print(f"epsilon_t = ({val_true} - {I})/ {val_true})*100")
        print(f"e_t is {epsilon_t}")
        
        
def simpson_1_3(f,a,b,n,val_true=False):
    i=0
    h = (b-a)/n
    print(f"h = ({b}-{a})/{n} = {h}")
    working = pd.DataFrame({"x":a,"f(x)":f(a),"state":-1},index=[i])
    for i in range(1,n+1):
        x = a+h*i
        working = working.append({"x":x,"f(x)":round(f(x),5),"state":i%2 if i != n else -1 },ignore_index=True)
    print(working)
    sum_even =np.sum(working[working.state == 0]['f(x)'])
    print(f"sum_even is {sum_even}")
    sum_odd =np.sum(working[working.state == 1]['f(x)'])
    print(f"sum_odd is {sum_odd}")
    
    I = (b-a)*((working['f(x)'][0]+4*sum_odd+2*sum_even+working['f(x)'][n])/(3*n))
    print(f"I = ({b}-{a})*(({working['f(x)'][0]}+4*{sum_odd}+2*{sum_even}+{working['f(x)'][n]})/(3*{n}))")
    print(round(I,5))
    
    if val_true:
        epsilon_t = ((val_true - I)/ val_true)*100
        epsilon_t = round(epsilon_t,5)
        print(f"epsilon_t = ({val_true} - {I})/ {val_true})*100")
        print(f"e_t is {epsilon_t}")
        
def simpson_3_8(f,a,b,n,val_true=False):
    i=0
    sigma=0
    h = (b-a)/n 
    print(f"h = ({b}-{a})/3 = {h}")
    working = pd.DataFrame({"x":a,"f(x)":f(a),"state":-1},index=[i])
    for i in range(1,n+1):
        x = a+h*i
        working = working.append({"x":x,"f(x)":round(f(x),5),"state":sigma if i!=n else -1},ignore_index=True)
        if i%3!=0:
            sigma +=1
        else:
            sigma=0
    print(working)

    sum_1 = np.sum(working[working.state == 0]['f(x)'])
    sum_2 = np.sum(working[working.state == 1]['f(x)'])
    sum_3 = np.sum(working[working.state == 2]['f(x)'])
    print(f"signma i =1.. --> {sum_1}")
    print(f"signma i =2.. --> {sum_2}")
    print(f"signma i =3... --> {sum_3}")
    I = ((3*h)/8)*((working['f(x)'][0]+3*sum_1+3*sum_2+2*sum_3+working['f(x)'][n]))
    print(f"I = ((3*{h})/8)*(({working['f(x)'][0]}+3*{sum_1}+3*{sum_2}+2*{sum_3}+{working['f(x)'][n]}))")
    print(I)
    
    if val_true:
        epsilon_t = ((val_true - I)/ val_true)*100
        epsilon_t = round(epsilon_t,5)
        print(f"epsilon_t = ({val_true} - {I})/ {val_true})*100")
        print(f"e_t is {epsilon_t}")

def mixed_simpson(f,a,b,n1,n3,val_true=False):
    i=0
    sigma=0
    h = (b-a)/(n1+n3) 
    print(f"h = ({b}-{a})/3 = {h}")
    working_1_3 = pd.DataFrame({"x":a,"f(x)":f(a),"state":-1},index=[i])
    
    for i in range(1,n1+1):
        x = a+h*i
        working_1_3 = working_1_3.append({"x":x,"f(x)":round(f(x),5),"state":i%2 if i != n1 else -1 },ignore_index=True)
    working_3_8 = pd.DataFrame({"x":working_1_3["x"][n1],"f(x)":working_1_3["f(x)"][n1],"state":-1},index=[i])
    for i in range(n1+1,n1+n3+1):
        x = a+h*i
        working_3_8 = working_3_8.append({"x":x,"f(x)":round(f(x),5),"state":sigma if i!=n1+n3 else -1},ignore_index=True)
        
        if sigma<2:
            sigma +=1
        else:
            sigma=0
        
    print("For Simpsons 1/3")
    print(working_1_3)
    sum_even =np.sum(working_1_3[working_1_3.state == 0]['f(x)'])
    print(f"sum_even is {sum_even}")
    sum_odd =np.sum(working_1_3[working_1_3.state == 1]['f(x)'])
    print(f"sum_odd is {sum_odd}")
    
    I_1_3 = h*((working_1_3['f(x)'][0]+4*sum_odd+2*sum_even+working_1_3['f(x)'][n1])/(3))
    print(f"I_1_3 = ({b}-{a})*(({working_1_3['f(x)'][0]}+4*{sum_odd}+2*{sum_even}+{working_1_3['f(x)'][n1]})/(3*{n1}))")
    print(round(I_1_3,5))
    
    print("For Simpson 3/8")
    print(working_3_8)
    sum_1 = np.sum(working_3_8[working_3_8.state == 0]['f(x)'])
    sum_2 = np.sum(working_3_8[working_3_8.state == 1]['f(x)'])
    sum_3 = np.sum(working_3_8[working_3_8.state == 2]['f(x)'])
    print(f"signma i =1.. --> {sum_1}")
    print(f"signma i =2.. --> {sum_2}")
    print(f"signma i =3... --> {sum_3}")
    I_3_8 = ((3*h)/8)*((working_3_8['f(x)'][0]+3*sum_1+3*sum_2+2*sum_3+working_3_8['f(x)'][n3]))
    print(f"I_3_8 = ((3*{h})/8)*(({working_3_8['f(x)'][0]}+3*{sum_1}+3*{sum_2}+2*{sum_3}+{working_3_8['f(x)'][n3]}))")
    print(I_3_8)
    
    print(f"\n\ntotal Area = {I_1_3} + {I_3_8} = {I_1_3 + I_3_8}")
    
    if val_true:
        epsilon_t = ((val_true - (I_1_3 + I_3_8))/ val_true)*100
        epsilon_t = round(epsilon_t,5)
        print(f"epsilon_t = ({val_true} - {I_1_3 + I_3_8})/ {val_true})*100")
        print(f"e_t is {epsilon_t}")