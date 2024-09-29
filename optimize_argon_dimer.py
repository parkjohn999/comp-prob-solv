from scipy.optimize import minimize
import matplotlib.pyplot as plt
import numpy as np
import os

def lennard_jones(r, epsilon=0.01, sigma=3.4):
    return 4*epsilon*(pow(sigma/r,12.0) - pow(sigma/r,6.0))

# Prevents the code from running when only the Lennard Jones function wants to be called
if __name__ == '__main__':
    # Minimizing the Lennard Jones Function
    minimization = minimize(
        fun=lennard_jones,
        x0=4,
        method="Nelder-Mead",
        tol=1e-6
    )

    print("Minimized Distance: ")
    print(minimization["x"][0])

    # Calculating the various Lennard Jones Potentials varying with distance
    r = np.arange(3,6,0.1).tolist()
    y = [lennard_jones(ret) for ret in r]

    # Writing into the xyz file
    file_dir = os.path.dirname(os.path.realpath('__file__'))
    file_name = os.path.join(file_dir, 'homework-2-1/dimer_geometry.xyz')
    File = open(file_name, 'w')
    toWrite = "2\nAr2\nAr 0 0 0\nAr " + str(round(minimization["x"][0])) + " 0 0" 
    File.write(toWrite)

    # Plotting the Distance v Voltage
    plt.plot(r, y)
    plt.axvline(x=minimization["x"])
    plt.xlabel("Distance (A)")
    plt.ylabel("Voltage (V)")
    plt.title("Lennard-Jones Potential Argon Dimer")
    plt.savefig('homework-2-1/ArgonDimerPotential.png')
    plt.show()

    