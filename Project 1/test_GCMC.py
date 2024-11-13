from GCMC import *
import numpy as np
import matplotlib.pyplot as plt

#print(compute_neighbor_indices(3))

# Parameters
size = 4
n_steps = 10000
mus_A = np.linspace(-0.2, 0, 7)
Ts = np.linspace(0.001, 0.019, 7)
params = []
for mu_A in mus_A:
    for T in Ts:
        params.append({
            'epsilon_A': -0.1,
            'epsilon_B': -0.1,
            'epsilon_AA': 0,
            'epsilon_BB': 0,
            'epsilon_AB': 0,
            'mu_A': mu_A,
            'mu_B': -0.1,
            'T': T  # Temperature (in units of k)
        })

#print(params[5])
sim = run_simulation(size,n_steps,params[5])

print(sim)

fig, ax = plt.subplots()

ax = plot_lattice(sim[0], ax, "the")

plt.show()

