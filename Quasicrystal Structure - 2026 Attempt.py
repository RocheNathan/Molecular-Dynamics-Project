#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 12:41:08 2026

@author: nathanroche
Title : Quasicrystal Project
"""

# Import the required packages
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

"SETUP THE INITIAL SYSTEM CONFIGURATION AND REQUIRED FUNCTIONS"
L = 18  # Container for the particles to sit in, a box of length L 
N = 150 # Number of particles inside of the box
T = 1.0 # Target temperature of the initial system

def initialise_positions(N , L):
# I have had to resort to using a starting lattice of particles to fit the box L
# I was getting a divide by 0 issue on r2 which was blowing up K    
    positions = []
    # How many particles to put in each row (sqrt of 150 is ~12.2, so 13)
    particles_per_side = int(np.ceil(np.sqrt(N)))
    
    # Calculate the gap between particles so they fit perfectly in the box
    spacing = L / particles_per_side
    
    for i in range(particles_per_side):
        for j in range(particles_per_side):
            if len(positions) < N:
                # Place x and y based on the current row (i) and column (j)
                # We add 'spacing/2' just to keep them away from the very edge
                x = i * spacing + (spacing / 2)
                y = j * spacing + (spacing / 2)
                positions.append([x, y])
                
    return np.array(positions)

# Assign the inital velocities based on Maxwell-Boltzmann distrobution
# Assumption here that particles of mass M = 1 and Boltzmann constant k = 1 
def initial_velocies(N , T):
    
    mu , sigma = 0.0 , np.sqrt(T)
    
    v_x = np.random.normal(mu, sigma, size = N)
    v_y = np.random.normal(mu , sigma, size = N)
    
    velocities = np.column_stack((v_x , v_y))
    
    # We will make the momentum of the system 0 so it does not drift off of the page/box
    velocities -= velocities.mean(axis = 0)
    
    return velocities

# Apply periodic boundary conditions for when the particles move if they 
# escape the box they return on the other side
# This function should be called after a time step
def apply_periodic_boundary_conditions(positions , L):
    
    # Using modulo will turns this flat box in to a infinate coordinate system on [0,L]
    positions %= L
    return positions

# To be defined but this is to calculate the impace of repulsion and attraction each particle i has on j
# Using reduced units we will take epsilon and sigma = 1.0
def calculate_forces(positions , L):
    total_potential_energy = 0.0
    epsilon = 1.0
    sigma = 1.0
    
    # 1. Initialize forces array as zeros
    N_local = len(positions)
    Forces = np.zeros((N_local,2))
    # 2. Loop through every pair (i, j)
    for i in range(0,N_local):
        for j in range(i+1 ,N_local):
    # 3. Calculate relative distance vector rij = pos[i] - pos[j]
                dx = positions[i,0] - positions[j,0]
                dy = positions[i,1] - positions[j,1]
                
    # 4. APPLY PERIODIC BOUNDARY CONDITIONS (The "Minimum Image Convention")
    #    rij = rij - L * np.round(rij / L)
                dx = dx - L * np.round(dx/L)
                dy = dy - L * np.round(dy/L)    
                
    # 5. Calculate scalar distance r = sqrt(dot(rij, rij))
    # Squared distance r2 is computationally better
                r2 = dx**2 + dy**2
                
                if r2 < 6.25:
    # 6. Calculate LJ Force magnitude
    # Computationally more efficient apparently...? 
                    r_inv6 = (1.0/r2)**(3.0)
                    r_inv12 = r_inv6**2
                
                    Force_magnitude = (48.0 * epsilon / r2) * (r_inv12 - 0.5*r_inv6)
    # 7. Add vector force to particle i, subtract from particle j
                    Forces[i,0] += Force_magnitude * dx
                    Forces[j,0] -= Force_magnitude * dx
                    Forces[i,1] += Force_magnitude * dy
                    Forces[j,1] -= Force_magnitude * dy
    
    # 7. Additionally calculating potential energy of the system - Sense Check
                    pair_potential = 4.0 * epsilon * (r_inv12 - r_inv6)
                    total_potential_energy += pair_potential
    
    return Forces, total_potential_energy

# Radial Distrobution Function g(r) is the mathematical confirmation of what shapes are appearing in the structure
def calculate_gr(positions, L , num_bins = 100):
    distances = []
    N = len(positions)
    max_r = L/2 # Can't calculate measure distances beyone half the box
    dr = max_r / num_bins
     
    for i in range(0,N):
        for j in range(i+1,N):
            dx = positions[i,0]-positions[j,0]
            dy = positions[i,1]-positions[j,1]
            
            dx = dx - L * np.round(dx/L)
            dy = dy - L * np.round(dy/L)  
            
            r = np.sqrt(dx**(2)+dy**(2))
            
            if r < max_r : distances.append(r)
    #Binning the distances together
    hist, bin_edges = np.histogram(distances, bins=num_bins, range=(0, max_r))
    
    # Calculate centers of bins for plotting
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # Normalization factor: density * area of ring * N
    rho = N / (L**2)
    # Area of ring = pi * ( (r+dr)^2 - r^2 )
    areas = np.pi * (bin_edges[1:]**2 - bin_edges[:-1]**2)
    
    # g(r) calculation 
    gr = hist / (areas * rho * N) * 2
    
    return bin_centers, gr
            


" SETUP THE SIMULATION ASPECT AND VELOCITY VERLET INTEGRATION "

# Simulation parameters
dt = 0.005 # This is the time step used, reduction in this will increase accuracy 
num_steps = 5000

trajectory = []         # To store positions for the animation
potential_history = []  # To store V
kinetic_history = []    # To store K
total_energy_history = []

#Annealing Parameters
inital_T = 1.0 # Starting temp - probably a liquid
final_T = 0.01 #Solid - this is for a super slow run

cooling_step = (inital_T-final_T)/num_steps #Proportianlly reduce the systems T to avoid flash freezing
current_target_T = inital_T

# 1. Get starting state of the system
positions = initialise_positions(N,L)
velocities = initial_velocies(N, T)

# 2. Get the initial force calculation completed
forces , potential_energy = calculate_forces(positions, L)

# 3. Time and Position advancement - Simulation 
for step in range(num_steps):
    # 3a. Half kick to velocity
    velocities = velocities +(0.5)*forces*dt
    
    # 3b. Positional drift with new velocity
    positions = positions +(velocities*dt)
    
    # 3c. Apply boundary conditions on new positions
    positions = apply_periodic_boundary_conditions(positions, L)
    
    # 3d. Force update
    forces , potential_energy = calculate_forces(positions, L)
    
    # 3e. Second half kick to velocity
    velocities = velocities +(0.5)*forces*dt
    
    if step % 10 == 0:
        trajectory.append(positions.copy()) # Taking snapshots of the position so we can plot
        
        # Calculate Kinetic Energy: K = 1/2 * m * v^2 (m=1 in reduced units)
        kinetic_energy = 0.5 * np.sum(velocities**2)
        
        # Save energies to lists
        potential_history.append(potential_energy)
        kinetic_history.append(kinetic_energy)
        total_energy_history.append(potential_energy + kinetic_energy)
    # 3f. Thermostat - Reducing the systems temperature to settle down in to a crystaline structure
        
    #What is the current temperature of the system. i.e velocity
    #Statistical mechanics aspect
        current_K = 0.5*np.sum(velocities**2)
        measured_T = current_K/N
        
        current_target_T -= (cooling_step*10)
        if current_target_T < final_T : current_target_T = final_T
        
        #Calculate rescaling factor
        # Ran in to division by 0 issue again, 1e^-10 should prevent this
        lamda_factor = np.sqrt(current_target_T/(current_target_T +(measured_T)*1**(-10)))
        
        velocities *= lamda_factor









"PLOTTING AND OUTPUTS"
plt.plot(total_energy_history, label="Total Energy")
plt.plot(potential_history, label="Potential")
plt.plot(kinetic_history, label="Kinetic")
plt.xlabel("Time Step")
plt.ylabel("Energy")
plt.title("Conservation of Energy in the system",)
plt.legend()
plt.show()


# Plot Frame 0 (Liquid/Grid)
plt.figure(figsize=(6,6))
plt.scatter(trajectory[0][:,0], trajectory[0][:,1])
plt.title("Initial State")
plt.show()

# Plot Final Frame 
plt.figure(figsize=(6,6))
plt.scatter(trajectory[-1][:,0], trajectory[-1][:,1])
plt.title("Final State (Annealed)")
plt.show()

plt.figure(figsize=(8,8))
final_positions = trajectory[-1]
plt.scatter(final_positions[:,0], final_positions[:,1], s=50, edgecolors='k', alpha=0.7)
plt.title("The 'Final Freeze' Structure")
plt.axis('equal')
plt.show()

# Plotting g(r) - Confirmation of shape
bins, gr_values = calculate_gr(trajectory[-1], L, num_bins=150)
plt.figure(figsize=(8, 5))

plt.plot(bins, gr_values, color='blue', linewidth=2)
plt.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5)
plt.xlabel(r'Distance $r$')
plt.ylabel(r'$g(r)$')
plt.title('Radial Distribution Function (Final State)')
plt.grid(alpha=0.3)
plt.show()

"""
Testing
pos = np.array([[10.0, 10.0], [10.5, 10.0]])
forces, energy = calculate_forces(pos, L)
print("Forces on Particle 0 and 1:")
print(forces)
print(f"Potential Energy: {energy}")

pos_pbc = np.array([[0.1, 10.0], [19.9, 10.0]])
L = 20.0
forces_pbc, energy_pbc = calculate_forces(pos_pbc, L)
print(forces_pbc)
print(energy_pbc)
"""
    



    
    
    
                
    
   