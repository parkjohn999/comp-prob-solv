import numpy as np

def initialize_lattice(size):
    return np.zeros((size,size))

def compute_neighbor_indices(size):
    neighbor_indices = []
    for x in range(0,size):
        for y in range(0,size):
            neighbors = [
                ((x-1) % size,y),
                ((x-1) % size,y),
                (x,(y-1) % size),
                (x,(y+1) % size)
            ]
            neighbor_indices[(x,y)] = neighbors
    return neighbor_indices

def calculate_interaction_energy(lattice, site, particle, neighbor_indices, epsilon_AA, epsilon_BB, epsilon_AB):
    x = site[0]
    y = site[1]
    interaction_energy = 0
    for neighbor in neighbor_indices:
        neighbor_particle = lattice[neighbor]
        if particle == 1:
            if neighbor_particle == 1:
                interaction_energy += (epsilon_AA)
            else:
                interaction_energy += (epsilon_AB)
        else:
            if neighbor_particle == 2:
                interaction_energy += epsilon_BB
            else:
                interaction_energy += epsilon_AB
    return interaction_energy