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

def attempt_move(lattice, N_A, N_B, N_empty, neighbor_indicies, params):
    size = len(lattice)
    N_sites = size * size
    beta = 1/params['T']
    epsilon_A, epsilon_B, epsilon_AA, epsilon_BB, epsilon_AB, mu_A, mu_B = params

    np.random.seed(42)
    if np.random.rand() > 0.5:
        # Add a particle
        if N_empty == 0:
            return (N_A, N_B, N_empty)
        # Random Empty Site in lattice
        x = np.random.randint(size - 1)
        y = np.random.randint(size - 1)
        while lattice[x,y] != 0:
            x = np.random.randint(size - 1)
            y = np.random.randint(size - 1)
        # Add Particle B
        particle = 2
        mu = mu_B
        epsilon = epsilon_B
        N_s = N_B
        if np.random.rand() > 0.5:
            # Add Particle A
            particle = 1
            mu = mu_A
            epsilon = epsilon_A
            N_s = N_A
        delta_E = epsilon + calculate_interaction_energy(lattice,(x,y),particle,neighbor_indicies,epsilon_AA,epsilon_BB,epsilon_AB)
        acc_prob = min[1, (N_empty) / (N_s + 1) * np.exp(-beta * (delta_E - mu))]
        r = np.random.rand()
        if r < acc_prob:
            lattice[(x,y)] = particle
            if particle == 1:
                N_A += 1
            else:
                N_B += 1
            N_empty -= 1
    else:
        # Remove a particle
        if N_sites - N_empty == 0:
            return (N_A, N_B, N_empty)
        # Random Occupied Site in lattice
        x = np.random.randint(size - 1)
        y = np.random.randint(size - 1)
        while not (lattice[x,y] == 0):
            x = np.random.randint(size - 1)
            y = np.random.randint(size - 1)
        particle = lattice[(x,y)]
        # Particle B
        mu = mu_B
        epsilon = epsilon_B
        N_s = N_B
        if particle == 1:
            mu = mu_A
            epsilon = epsilon_A
            N_s = N_A
        delta_E = -epsilon - calculate_interaction_energy(lattice,(x,y),particle,neighbor_indicies,epsilon_AA,epsilon_BB,epsilon_AB)
        acc_prob = min[1, N_s / (N_empty + 1) * np.exp(-beta * (delta_E + mu))]
        r = np.random.rand()
        if r < acc_prob:
            lattice[(x,y)] = 0
            if particle == 1:
                N_A -= 1
            else:
                N_B -= 1
            N_empty += 1
    return (N_A, N_B, N_empty)

def run_simulation(size, n_steps, params):
    lattice = initialize_lattice(size)
    neighbor_indices = compute_neighbor_indices(size)
    N_sites = size * size
    N_A = 0
    N_B = 0
    N_empty = N_sites
    
    coverage_A = np.zeros(n_steps)
    coverage_B = np.zeros(n_steps)

    for step in range(n_steps):
        N_A, N_B, N_empty = attempt_move(lattice, N_A, N_B, N_empty, neighbor_indices, params)
        coverage_A[step] = N_A / N_sites
        coverage_B[step] = N_B / N_sites
    
    return (lattice, coverage_A, coverage_B)

def plot_lattice(lattice, ax, title):
    size = len(lattice)
    for x in range(size):
        for y in range(0,size):
            if lattice[x,y] == 1:
                ax.plot(x + 0.5, y + 0.5, c="red")
            elif lattice[x,y] == 2:
                ax.plot(x + 0.5, y + 0.5, c="blue")
    ax.set_xlim([0,size])
    ax.set_ylim([0,size])
    ax.xticks([])
    ax.yticks([])
    ax.minorticks_on()
    ax.title(title)
    return ax