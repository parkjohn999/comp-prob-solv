import scipy
import numpy as np
import scipy.integrate

R = 8.314

def to_int(V, p_const, aindex):
    return p_const/(pow(V, aindex))

def compute_work_adiabatic(Vi, Vf, aindex=1.4, n=1, T=300):
    # Constant Computation
    # Initial Pressure Computation using PV = nRT
    Pi = (n*R*T)/Vi
    # Solve for constant
    p_const = Pi * pow(Vi, aindex)
    
    # Integrate
    x = np.linspace(Vi, Vf, num=100)
    y = to_int(x, p_const, aindex)
    return scipy.integrate.trapezoid(y,x=x)