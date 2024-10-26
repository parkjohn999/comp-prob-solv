import math
import numpy as np
from scipy.stats import expon

np.random.seed(42)

def polar_to_cartesian (r, theta, alpha=0):
    """
    Converts from Polar Coordinates into Cartesian Coordinates

    r: radius
    theta: horizontal angle from the x-axis in radians
    alpha: vertical angle from the z-axis in radians

    alpha defaults to 0 since in our assumption the 2p orbital is aligned with the Z axis
    
    Returns tuple (x,y,z)
    """
    return (
        r*np.sin(theta)*np.cos(alpha),
        r*np.sin(theta)*np.sin(alpha),
        r*np.cos(theta)
        )

def cartesian_to_polar(x, y, z):
    """
    Converts from Cartesian Coordinates into Polar Coordinates
    
    Returns tuple (r, theta, alpha) in radians
    """
    r=np.sqrt(np.power(x,2)+np.power(y,2)+np.power(z,2))
    return (
        r,
        np.arcsin(z/r),
        np.arctan(y/x)
        )

def psi_2p_z(x,y,z, alpha_0 = 1):
    """
    Computes the hydrogen atom's 2p orbital function
    Assumes 2p orbital is oriented on z-axis

    x: distance along x-axis
    y: distance along y-axis
    z: distance along z-axis
    alpha_0: Bohr Radius
    """
    # Convert into Polar Coordinates
    polar_coords = cartesian_to_polar(x,y,z)
    r = polar_coords[0]
    theta = polar_coords[1]

    # Calculation
    return 1/(4*np.sqrt(2*np.pi)*np.power(alpha_0, 1.5))*(r/alpha_0)*np.cos(theta)*np.exp(-r/(2*alpha_0))

# Task 1.1.2

# Constants
R = 2
L = 20
# My Computer is old so I couldn't get 10^8 to run without crashing
N_POINTS = [int(math.pow(10, i)) for i in range(2,8)]

def random_sampling_integral():
    """
    Performs the Integral using random sampling

    Returns S(R)
    """

    # Collection Array
    S_Rs = []

    # Loop over so we generate the necessary amount of points
    for n in N_POINTS:
        # Random Point Generation with L
        x = np.random.uniform(-L,L,n)
        y = np.random.uniform(-L,L,n)
        z = np.random.uniform(-L,L,n)
        # Integrand Calculation
        integrand = psi_2p_z(x,y,z + R/2)*psi_2p_z(x,y,z-R/2)
        
        # S(R) Calculation and append to collection array
        S_R = pow(2*L, 3) * np.average(integrand)
        S_Rs.append(S_R)

    return S_Rs

# Task 1.1.3

def importance_sampling_integral():
    """
    Performs the integral using the exponential distribution
    
    Returns S(R)
    """

    # Collection Array
    S_Rs = []

    # Loop over so we generate the necessary amount of points
    for n in N_POINTS:
        # Point Generation
        x = expon.rvs(size=n,loc=2,scale=1/(0.0053))
        y = expon.rvs(size=n,loc=2,scale=1/(0.0053))
        z = expon.rvs(size=n,loc=2,scale=1/(0.0053))
        # Integrand Calculation
        numer = psi_2p_z(x,y,z + R/2)*psi_2p_z(x,y,z-R/2)
        denom = expon.pdf(x,loc=2,scale=1/(0.0053)) * expon.pdf(y,loc=2,scale=1/(0.0053)) * expon.pdf(z,loc=2,scale=1/(0.0053))
        integrand = numer / denom

        # S(R) Calculation and append to collection array
        S_R = np.average(integrand)
        S_Rs.append(S_R)

    return S_Rs

# Task 1.1.4
def importance_sampling_integral_R():
    """
    Performs the integral using the exponential distribution
    
    Returns S(R)
    """

    # Collection Array
    S_Rs = []
    n = pow(10,4)
    Rs = np.linspace(0.5, 20, 40)

    # Loop over so we generate the necessary amount of points
    for r in Rs:
        # Point Generation
        x = expon.rvs(size=n,loc=2,scale=1/(0.0053))
        y = expon.rvs(size=n,loc=2,scale=1/(0.0053))
        z = expon.rvs(size=n,loc=2,scale=1/(0.0053))
        # Integrand Calculation
        numer = psi_2p_z(x,y,z + r/2)*psi_2p_z(x,y,z-r/2)
        denom = expon.pdf(x,loc=2,scale=1/(0.0053)) * expon.pdf(y,loc=2,scale=1/(0.0053)) * expon.pdf(z,loc=2,scale=1/(0.0053))
        integrand = numer / denom

        # S(R) Calculation and append to collection array
        S_R = np.average(integrand)
        S_Rs.append(S_R)

    return S_Rs