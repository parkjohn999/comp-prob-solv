from lecture7_functions import *
from lecture8_functions import *
import pandas as pd
import matplotlib.pyplot as plt

# Read in the csv file
df = pd.read_csv("trouton.csv")
# Convert Hv into J/mol
Hv = df["H_v (kcal/mol)"] * 4184
TB = df["T_B (K)"]
classes = df["Class"]

# Compute the slope and intercept for the linear regression model
my_ols = ols(TB, Hv)
slope = my_ols[0]
intercept = my_ols[1]

# Print the results
print("SLOPE:")
print(slope)
print("INTERCEPT: ")
print(intercept)

# Interpretation
# The slope here is the entropy of vaporization 
# It documents for an incease in 1 K in boiling point, what hte corresponding change is the enthalpy of vaporization

# The slope is 103.85 J/mol K so Trouton's Rule works reasonably well
# This is only a 15% difference

# 95% Confidence Interval calculation
my_line = slope * TB + intercept
residuals = Hv - my_line
# Slope
slope_confidence = confidence_interval_slope(TB,residuals,0.95)
print("SLOPE CONFIDENCE")
print(slope_confidence)
# Intercept
intercept_confidence = confidence_interval_intercept(TB, residuals, 0.95)
print("INTERCEPT CONFIDENCE")
print(intercept_confidence)

# Preperation to color points
my_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22']
possible_classes = ["Perfect liquids", "Imperfect liquids", "Liquids subject to quantum effects", "Metals"]
plot_colors = []
for c in classes:
    plot_colors.append(my_colors[possible_classes.index(c)])

# Plot H_v vs T_B with the linear regression line
plt.scatter(TB, Hv, c=plot_colors)
plt.plot(TB, my_line)
equation = f'H_v = {slope:.2f} ± {slope_confidence:.2f} * T_B + {intercept:.2f} ± {intercept_confidence:.2f}'
plt.title("Trouton's Rule\n" + equation)
plt.xlabel("T_B (K)")
plt.ylabel("H_v (J/(mol K))")
plt.savefig("homework-3-1/TroutonsRule.png")
plt.show()
