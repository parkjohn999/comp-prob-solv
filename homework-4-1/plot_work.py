from compute_work_adiabatic import *
from compute_work_isothermal import *
import scipy
import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd

# Constants
Vi = 0.1
Vf = 3.0*Vi

# Setup Lists for Plotting
x = np.linspace(Vi, Vf, num=100)
adiabatics = [compute_work_adiabatic(Vi, v) for v in x]
isothermals = [compute_work_isothermal(Vi, v) for v in x]

# Write to csv
df = pd.DataFrame()
df['Final Volume'] = x
df['Adiabatic Work'] = adiabatics
df['Isothermal Work'] = isothermals
df.to_csv('volumeAndWork.csv', index=False)

# Plotting
plt.plot(x,adiabatics,label="Adiabatic Work")
plt.plot(x,isothermals,label="Isothermal Work")
plt.xlabel("Final Volume (m^3)")
plt.ylabel("Work (J)")
plt.legend(loc="upper left")
plt.title("Isothermal & Adiabatic Work")
plt.show()
