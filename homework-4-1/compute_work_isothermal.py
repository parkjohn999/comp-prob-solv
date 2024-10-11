import scipy
import numpy as np

R = 8.314

def to_int(V, n, T):
    return (n*R*T)/V

def compute_work_isothermal(Vi, Vf, n=1, T=300):
    x = np.linspace(Vi, Vf, num=100)
    y = to_int(x, n, T)
    return scipy.integrate.trapezoid(y,x=x)