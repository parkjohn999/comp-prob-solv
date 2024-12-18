from analysis_functions import *
import matplotlib.pyplot as plt
from MD_simulation import *

# Simulation parameters
dt = 0.01  # Time step
total_steps = 10000  # Number of steps
box_size = 100.0  # Size of the cubic box
k = 1.0  # Spring constant
mass = 1.0  # Particle mass
r0 = 1.0  # Equilibrium bond length
target_temperature = 0.1  # Target temperature
rescale_interval = 100  # Steps between velocity rescaling
n_particles = 25  # Number of particles
epsilon_repulsive = 1.0  # Depth of repulsive LJ potential
epsilon_attractive = 0.5  # Depth of attractive LJ potential
sigma = 1.0  # LJ potential parameter

# Arrays to store properties
temperatures = np.linspace(0.1, 1.0, 10)
Rg_values = []
Ree_values = []
potential_energies = []

for T in temperatures:
    # Set target temperature
    target_temperature = T
    # (Re-initialize positions and velocities)
    positions = initialize_chain(n_particles, box_size, r0)
    velocities = initialize_velocities(n_particles, target_temperature, mass)
    potential_energy_array = np.zeros_like(positions)
    # (Run simulation)

    for step in range(total_steps):
        forces_harmonic = compute_harmonic_forces(positions, k, r0, box_size)
        forces_lennard_jones = compute_lennard_jones_forces(positions, epsilon_repulsive, epsilon_attractive, sigma, box_size)
        total_forces = forces_harmonic + forces_lennard_jones
        potential_energy_array = total_forces
        
        # Integrate equations of motion
        positions, velocities, total_forces = velocity_verlet(positions, velocities, total_forces, dt, mass, epsilon_repulsive, epsilon_attractive, sigma, box_size, k, r0)
        
        # Apply thermostat
        if step % rescale_interval == 0:
            velocities = rescale_velocities(velocities, target_temperature, mass, n_particles)

    # Compute properties
    Rg = calculate_radius_of_gyration(positions)
    Ree = calculate_end_to_end_distance(positions)
    Rg_values.append(Rg)
    Ree_values.append(Ree)
    potential_energies.append(np.mean(potential_energy_array))

# Plotting
plt.figure()
plt.plot(temperatures, Rg_values, label='Radius of Gyration')
plt.xlabel('Temperature')
plt.ylabel('Radius of Gyration')
plt.title('Radius of Gyration vs Temperature')
plt.legend()
plt.show()

plt.figure()
plt.plot(temperatures, Ree_values, label='End-to-End Distance')
plt.xlabel('Temperature')
plt.ylabel('End-to-End Distance')
plt.title('End-to-End Distance vs Temperature')
plt.legend()
plt.show()

plt.figure()
plt.plot(temperatures, potential_energies, label='Potential Energy')
plt.xlabel('Temperature')
plt.ylabel('Potential Energy')
plt.title('Potential Energy vs Temperature')
plt.legend()
plt.show()