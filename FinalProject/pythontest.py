import numpy as np 
from numdifftools import Hessian

def calculate_LJ (rij, sigma, epsilon):
    r12 = (sigma/rij)**12
    r6 = (sigma/rij)**6
    return epsilon*(r12 - (2*r6)) 

def calculate_energy(coords, sigma, epsilon):
    num_atoms = 4 
    energy = 0.0 
    for i in range(num_atoms):
        for j in range(i+1, num_atoms):
            rij = np.linalg.norm(coords[i] - coords[j])
            energy += calculate_LJ(rij, sigma[i,j], epsilon[i,j])
    return energy

# Test cases 
num_atoms = 4
test_coords = np.array([[4.4768, 0, 0], [ 1.7582, 0, 0], [-1.7582, 0, 0], [-4.4768, 0, 0]])
test_sigma = np.array([[2.55, 2.72, 2.72, 2.55], [2.72, 3.52, 3.52, 2.73], [2.72, 3.52, 3.52, 2.72], [2.55, 2.72, 2.72, 2.55]])
test_epsilon = np.array([[0.015, 0.05, 0.05, 0.015], [0.05, 0.066, 0.066, 0.05], [0.05, 0.066, 0.066, 0.05], [0.015, 0.05, 0.05, 0.015]])


def energy_wrapper(coords):
    return calculate_energy(coords.reshape((num_atoms, 3)), test_sigma, test_epsilon)

hessian_calculator = Hessian(energy_wrapper)
hessian = hessian_calculator(test_coords.flatten())

print(hessian)

