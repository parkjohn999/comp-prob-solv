# Import necessary libraries
import numpy as np
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
n_particles = 20  # Number of particles
epsilon_repulsive = 1.0  # Depth of repulsive LJ potential
epsilon_attractive = 0.5  # Depth of attractive LJ potential
sigma = 1.0  # LJ potential parameter

# Initialize positions and velocities
positions = initialize_chain(n_particles, box_size, r0)
velocities = initialize_velocities(n_particles, target_temperature, mass)

forces_array = []

# Simulation loop
for step in range(total_steps):
    # Compute forces
    forces_harmonic = compute_harmonic_forces(positions, k, r0, box_size)
    forces_lennard_jones = compute_lennard_jones_forces(positions, epsilon_repulsive, epsilon_attractive, sigma, box_size)
    total_forces = forces_harmonic + forces_lennard_jones
    
    # Integrate equations of motion
    positions, velocities, total_forces = velocity_verlet(positions, velocities, total_forces, dt, mass, epsilon_repulsive, epsilon_attractive, sigma, box_size, k, r0)
    
    # Apply thermostat
    if step % rescale_interval == 0:
        velocities = rescale_velocities(velocities, target_temperature, mass, n_particles)

    forces_array.append(np.mean(total_forces))

#plt.plot(range(total_steps), forces_array)
#plt.show()
plt.plot(range(total_steps), forces_array)
plt.show()