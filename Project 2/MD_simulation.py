import numpy as np
import math

# Constants
k_B = 1.380649e-23

def random_unit_vector():
    vector = np.random.normal(0, 1, 3)
    magnitude = np.linalg.norm(vector)
    unit_vector = vector / magnitude
    return unit_vector

def apply_pbc(position, box_size):
    return position % box_size

def minimum_image(displacement, box_size):
    return displacement - box_size * np.round(displacement / box_size)

def initialize_chain(n_particles, box_size, r0):
    positions = np.zeros((n_particles, 3))
    current_position = [box_size/2, box_size/2, box_size/2]
    positions[0] = current_position
    for i in range(1, n_particles - 1):
        direction = random_unit_vector()
        next_position = current_position + r0 * direction
        positions[i] = apply_pbc(next_position, box_size)
        current_position = positions[i]
    return positions

def initialize_velocities(n_particles, target_temperature, mass):
    velocities = np.random.normal(0, np.sqrt(k_B * target_temperature / mass), (n_particles, 3))
    velocities -= np.mean(velocities)
    return velocities

def compute_harmonic_forces(positions, k, r0, box_size):
    forces = np.zeros_like(positions)
    for i in range(0, len(positions) - 1):
        displacement = positions[i+1] - positions[i]
        displacement = minimum_image(displacement, box_size)
        distance = np.linalg.norm(displacement)
        force_magnitude = -k * (distance - r0)
        force = force_magnitude * (displacement / distance)
        forces[i] -= force
        forces[i + 1] += force

    return forces

def compute_lennard_jones_forces(positions, epsilon_repulsive, epsilon_attractive, sigma, box_size):
    forces = np.zeros_like(positions)
    for i in range(0,len(positions)):
        for j in range(i+2, len(positions)):
            # The idea is to model beads seperated by one spacer with repulsive LJ
            # and model beads seperated by more than one spacer with attractive LJ
            displacement = positions[j] - positions[i]
            displacement = minimum_image(displacement, box_size)
            distance = np.linalg.norm(displacement)
            if abs(i-j) == 2:
                # Repulsive Senerio when seperated by one spacer
                cutoff = 1.122 * sigma
                if distance < cutoff:
                    force_magnitude = 4 * epsilon_repulsive * (np.power((sigma / distance),12)-np.power(sigma/distance,6)+0.25)
                    force = force_magnitude * (displacement / distance)
                    forces[i] -= force
                    forces[j] += force
            else:
                # Attractive Senerio when seperated by more than one spacer
                force_magnitude = 4 * epsilon_attractive * (np.power((sigma / distance),12)-np.power(sigma/distance,6))
                force = force_magnitude * (displacement / distance)
                forces[i] -= force
                forces[j] += force

    return forces    

def compute_forces(positions, epsilon_repulsive, epsilon_attractive, sigma, box_size, k, r0):
    # Compute forces
    forces_harmonic = compute_harmonic_forces(positions, k, r0, box_size)
    forces_lennard_jones = compute_lennard_jones_forces(positions, epsilon_repulsive, epsilon_attractive, sigma, box_size)
    return forces_harmonic + forces_lennard_jones

def velocity_verlet(positions, velocities, forces, dt, mass, epsilon_repulsive, epsilon_attractive, sigma, box_size, k, r0):
    velocities += 0.5 * forces / mass * dt
    positions += velocities * dt
    positions = apply_pbc(positions, box_size)
    forces_new = compute_forces(positions, epsilon_repulsive, epsilon_attractive, sigma, box_size, k, r0)
    velocities += 0.5 * forces_new / mass * dt
    return positions, velocities, forces_new

def rescale_velocities(velocities, target_temperature, mass, n_particles):
    kinetic_energy = 0.5 * mass * np.sum(np.power(np.linalg.norm(velocities, axis=1), 2))
    current_temperature = (2/3) * kinetic_energy / (n_particles * k_B)
    scaling_factor = np.sqrt(target_temperature / current_temperature)
    velocities *= scaling_factor
    return velocities
