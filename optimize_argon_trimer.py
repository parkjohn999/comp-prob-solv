from scipy.optimize import minimize
import matplotlib.pyplot as plt
import numpy as np
from optimize_argon_dimer import lennard_jones
import math
import os

def compute_bond_length(coord1, coord2):
    """
    Computes the bond lengths for a molecule given a set of coordinates.

    Parameters:
    coord1 (list): coordinates of a molecule.
    coord2 (list): coordinates of a molecule.

    Returns:
    float: The bond length.
    """
    # Check coord1 and coord2 are the same size
    if len(coord1) != len(coord2):
        print("ERROR in compute_bond_length coordinates")
    squaredSum = 0
    for i in range(len(coord1)):
        squaredSum += pow((coord2[i] - coord1[i]),2.0)
    sqrtSum = math.sqrt(squaredSum)
    # Check bond length is not too long
    if sqrtSum > 2:
        print("ERROR in compute_bond_length: bond length is too long")
    return sqrtSum

def vectorMaker(coord1, coord2):
    # Check coord1 and coord2 are the same size
    if len(coord1) != len(coord2):
        print("ERROR in vectorMaker coordinates")
    vector = np.zeros(len(coord1))
    for i in range(len(vector)):
        vector[i] = coord2[i] - coord1[i]
    return vector

def dotProd(vec1, vec2):
    # Check vec1 and vec2 are the same size
    if len(vec1) != len(vec2):
        print("ERROR in dotProd vectors")
    vecSum = 0
    for i in range(len(vec1)):
        vecSum += vec1[i] * vec2[i]
    return vecSum

def vecLength(vec1):
    vecSum = 0
    for i in range(len(vec1)):
        vecSum += pow(vec1[i],2.0)
    return math.sqrt(vecSum)

def compute_bond_angle(coord1, coord2, coord3):
    """
    Computes the bond angles for a molecule given a set of coordinates.

    Parameters:
    coord1 (list): coordinates of a molecule.
    coord2 (list): coordinates of a molecule.
    coord3 (list): coordinates of a molecule.

    Returns:
    float: The bond angle.
    """
    AB = vectorMaker(coord1, coord2)
    BC = vectorMaker(coord2, coord3)
    return math.degrees(math.acos((dotProd(AB,BC))/(vecLength(AB)*vecLength(BC))))

def sum_lennard_jones(params):
    r12, x3, y3 = params
    r13 = compute_bond_length([0,0,0],[x3, y3, 0])
    r23 = compute_bond_length([r12, 0, 0], [x3, y3, 0])
    return lennard_jones(r13) + lennard_jones(r12) + lennard_jones(r23)

# Minimizing the new Lennard Jones Function
result = minimize(sum_lennard_jones, [1,1,1])

print("Minimized Distances")
x3 = result["x"][1]
y3 = result["x"][2]
print("r12")
print(result["x"][0])
print("r13")
print(compute_bond_length([0,0,0],[x3, y3, 0]))
print("r23")
print(compute_bond_length([result["x"][0], 0, 0], [x3, y3, 0]))

# The Shape appears to be an Isosceles Triangle

print("Minimized Angles")
print("Atom 1 angle:")
print(compute_bond_angle([0,0,0], [x3, y3, 0], [result["x"][0], 0, 0]))
print("Atom 2 angle")
print(compute_bond_angle([result["x"][0], 0, 0], [0,0,0], [x3, y3, 0]))
print("Atom 3 angle")
print(compute_bond_angle([x3, y3, 0], [result["x"][0], 0, 0], [0,0,0]))

# Writing to the XYZ File
file_dir = os.path.dirname(os.path.realpath('__file__'))
file_name = os.path.join(file_dir, 'homework-2-1/trimer_geometry.xyz')
File = open(file_name, 'w')
toWrite = "3\nAr3\nAr 0 0 0\nAr " + str(round(result["x"][0])) + " 0 0\nAr " + str(round(x3)) + " " + str(round(y3)) + " 0"
File.write(toWrite)