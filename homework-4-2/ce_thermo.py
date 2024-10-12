import scipy
import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants
import pandas as pd

# My Constants
Ti = 300
Tf = 2000
T =np.linspace(Ti,Tf,100)
k_B = scipy.constants.k/scipy.constants.eV

# Repeatable functions so I don't have to write the same code 3 times
def internal_energy_calculation(Z, Temp=T):
    return -1 * np.gradient(np.log(Z), 1 / (k_B * Temp))

def free_energy_calculation(Z, Temp = T):
    return -1 * k_B*Temp*np.log(Z)

def entropy_calculation(F, Temp=T):
    return -1 * np.gradient(F,Temp)

# Functions for tasks
def isolated_partition(Temp=T):
    # Since the energies are zero, the partition function just sums up 1
    return [1*14 for t in Temp]

def soc_partition(Temp=T):
    # 6 zero energy states and 8 0.28 energy states
    return [6 + 8*(math.exp(-0.28/(k_B*t))) for t in Temp]

def cfs_partition(Temp=T):
    # 4 zero energy, 2 0.12 energy, 2 0.25 energy, 4 0.32 energy, 2 0.46 energy
    return [(4.0 + 2.0*math.exp(-0.12/(k_B*t)) + 2.0*math.exp(-0.25/(k_B*t)) + 4.0*math.exp(-0.32/(k_B*t)) + 2.0*math.exp(-0.46/(k_B*t))) for t in Temp]

# Printing Results for Isolated Ce(3)
isolated_Z = isolated_partition()
print("Isolated Internal Energy")
print(internal_energy_calculation(isolated_Z)[0])
print("Isolated Free Energy")
isolated_free_energy = free_energy_calculation(isolated_Z)
print(isolated_free_energy[0])
print("Isolated Entropy")
print(entropy_calculation(isolated_free_energy)[0])
print("\n")

# Printing Results for Ce(3) with SOC
SOC_Z = soc_partition()
print("SOC Internal Energy")
print(internal_energy_calculation(SOC_Z)[3])
print("SOC Free Energy")
SOC_free_energy = free_energy_calculation(SOC_Z)
print(SOC_free_energy[0])
print("SOC Entropy")
print(entropy_calculation(SOC_free_energy)[0])
print("\n")

# Printing Results for Ce(3) with SOC and CFS
cfs_Z = cfs_partition()
print("cfs Internal Energy")
print(internal_energy_calculation(cfs_Z)[0])
print("cfs Free Energy")
cfs_free_energy = free_energy_calculation(cfs_Z)
print(cfs_free_energy[0])
print("cfs Entropy")
print(entropy_calculation(cfs_free_energy)[0])

# CSV file
df = pd.DataFrame()
df["Isolated Ce(3) Internal Energy"] = internal_energy_calculation(isolated_Z)
df["Isolated Ce(3) Free Energy"] = free_energy_calculation(isolated_Z)
df["Isolated Ce(3) Entropy"] = entropy_calculation(isolated_free_energy)
df["Ce(3) with SOC Internal Energy"] = internal_energy_calculation(SOC_Z)
df["Ce(3) with SOC Free Energy"] = free_energy_calculation(SOC_Z)
df["Ce(3) with SOC Entropy"] = entropy_calculation(SOC_free_energy)
df["Ce(3) with SOC and CFS Internal Energy"] = internal_energy_calculation(cfs_Z)
df["Ce(3) with SOC and CFS Free Energy"] = free_energy_calculation(cfs_Z)
df["Ce(3) with SOC and CFS Entropy"] = entropy_calculation(cfs_free_energy)
df.to_csv("ThermoProperties.csv",index=False)

# Plotting
fig, (ax1, ax2, ax3) = plt.subplots(1,3,figsize=(12,5))

# Plot Internal Energy
ax1.plot(T, internal_energy_calculation(isolated_Z), label="Isolated Ce(3)")
ax1.plot(T, internal_energy_calculation(SOC_Z), label="Ce(3) with SOC")
ax1.plot(T, internal_energy_calculation(cfs_Z), label="Ce(3) with SOC and CFS")
ax1.set_xlabel("Temperature (K)")
ax1.set_ylabel("Internal Energy (eV)")
ax1.set_title("Internal Energy")
ax1.legend()

# Plot Free Energy
ax2.plot(T, free_energy_calculation(isolated_Z), label="Isolated Ce(3)")
ax2.plot(T, free_energy_calculation(SOC_Z), label="Ce(3) with SOC")
ax2.plot(T, free_energy_calculation(cfs_Z), label="Ce(3) with SOC and CFS")
ax2.set_xlabel("Temperature (K)")
ax2.set_ylabel("Free Energy (eV)")
ax2.set_title("Free Energy")
ax2.legend()

# Plot Entropy
ax3.plot(T, entropy_calculation(isolated_Z), label="Isolated Ce(3)")
ax3.plot(T, entropy_calculation(SOC_Z), label="Ce(3) with SOC")
ax3.plot(T, entropy_calculation(cfs_Z), label="Ce(3) with SOC and CFS")
ax3.set_xlabel("Temperature (K)")
ax3.set_ylabel("Entropy (eV)")
ax3.set_title("Entropy")
ax3.legend()

plt.savefig("ce_thermo.png")
plt.show()
