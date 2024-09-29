from scipy.optimize import minimize
from scipy.integrate import trapezoid
import matplotlib.pyplot as plt
import numpy as np
from optimize_argon_dimer import lennard_jones
import math
import os
import csv

def hard_sphere_potential(r, sigma=3.4):
    if r < sigma:
        return 1000
    
    return 0

def square_well_potential(r, sigma=3.4, labda=1.5, epsilon=0.01):
    if r < sigma:
        return 1000
    
    if r >= (labda * sigma):
        return 0
    
    return epsilon

# Constants
BOLTZMAN_CONST = 1.380 * pow(10,(-23))
PI = 3.14
AVOGADRO_NUMBER = 6.0221 * pow(10,23)

def to_int(r, temperature, potential_func):
    try:
        return (pow(r,2.0))*(math.exp(-potential_func(r)/(BOLTZMAN_CONST*temperature))-1)
    except:
        return 0

def B2V(temperature, potential_func, sigma=3.4, labda=1.5, epsilon=0.01):
    r = np.arange(0.01, 5*sigma, 0.01).tolist()
    y = [to_int(res, temperature, potential_func) for res in r]
    inte = trapezoid(y, r)
    return -2*PI*AVOGADRO_NUMBER*inte

# Print the B2V values at 100K
print("HARD SPHERE B2V")
print(B2V(100, hard_sphere_potential))

print("SQUARE WELL B2V")
print(B2V(100, square_well_potential))

print("LENNARD-JONES B2V")
print(B2V(100, lennard_jones))

# Calculate the necessary B2V constants varying with temperature
temps = range(100, 800)
B2V_hard_values = [B2V(temp, hard_sphere_potential) for temp in temps]
B2V_square_values = [B2V(temp, square_well_potential) for temp in temps]
B2V_lennard_values = [B2V(temp, lennard_jones) for temp in temps]

# Write the values into a csv file
with open("homework-2-2/B2V.csv", "w") as f:
    f.write("Temperature,Hard Sphere,SQUARE WELL,LENNARD-JONES\n")
    for i in range(len(temps)):
        f.write(str(temps[i]) + "," + str(B2V_hard_values[i]) + "," + str(B2V_square_values[i]) + "," + str(B2V_square_values[i]))

# Plot the values and save the resulting plot
plt.plot(temps, B2V_hard_values)
plt.plot(temps, B2V_square_values)
plt.plot(temps, B2V_lennard_values)
plt.xlabel("Temperature (K)")
plt.ylabel("B2V")
plt.title("B2V Values Versus Temperature")
plt.savefig('homework-2-2/B2V.png')
plt.show()