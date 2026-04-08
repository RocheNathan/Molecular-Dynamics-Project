# Molecular-Dynamics-Project

# Quasicrystal Project: Phase 1 (Lennard-Jones Annealing)

This project implements a 2D Molecular Dynamics (MD) simulation engine in Python. The current phase focuses on simulating the transition of a particle system from a high-energy gas/liquid state into a crystalline hexagonal lattice using a **Lennard-Jones (12-6) Potential** and simulated annealing.

## 🚀 Overview

The simulation uses the **Velocity Verlet Integration** method to solve the equations of motion for $N$ particles in a square box with **Periodic Boundary Conditions (PBC)**. To achieve a stable crystalline structure, the system undergoes a controlled cooling process (thermostat rescaling).

### Key Features:
* **Velocity Verlet Integration:** A second-order integrator to ensure energy stability.
* **Minimum Image Convention:** Accurate force calculations across periodic boundaries.
* **Simulated Annealing:** A velocity-rescaling thermostat that gradually lowers the system temperature.
* **Structural Analysis:** Implementation of the **Radial Distribution Function $g(r)$** to mathematically confirm the formation of hexagonal symmetry.

---

## 🛠 Simulation Parameters (Current)

| Parameter | Value | Description |
| :--- | :--- | :--- |
| `N` | 150 | Number of particles |
| `L` | 18 | Box side length |
| `dt` | 0.005 | Integration time step (reduced units) |
| `initial_T` | 1.0 | Starting system temperature |
| `final_T` | 0.01 | Target "freeze" temperature |
| `num_steps` | 5000 | Total simulation steps |

---

## 📊 Structural Analysis: $g(r)$

The **Radial Distribution Function** is used to identify the phase of the system. In the current hexagonal crystal state, we expect to see distinct peaks at specific distance ratios:
* **$1.0 \times a$**: Nearest neighbors.
* **$\sqrt{3} \times a$**: Second neighbors (hexagonal signature).
* **$2.0 \times a$**: Third neighbors.



---

## 📁 Project Structure

* `initialise_positions`: Sets up a stable starting lattice to prevent initial overlaps.
* `calculate_forces`: Computes the 12-6 Lennard-Jones forces and potential energy.
* `apply_periodic_boundary_conditions`: Handles particle wrapping at box edges.
* `calculate_gr`: Post-processing function to generate the $g(r)$ histogram.

---

## 🛠 Installation & Usage

1.  **Requirements:**
    * Python 3.x
    * NumPy
    * Matplotlib

2.  **Running the simulation:**
    ```bash
    python quasicrystal_project.py
    ```

3.  **Outputs:**
    * Energy Conservation Plot (Kinetic, Potential, Total).
    * Initial vs. Final state scatter plots.
    * Radial Distribution Function ($g(r)$) plot confirming the lattice type.

---

## 🚧 Future Work: The LJG Branch

The next phase of the project involves branching from this stable MD engine to implement a **Lennard-Jones-Gauss (LJG)** potential. By introducing a secondary attraction well at an "incommensurate" distance (e.g., $1.5 \times \sigma$), the system will be frustrated, preventing standard hexagonal packing and encouraging the formation of **quasicrystalline** symmetries (5-fold or 10-fold).



---
**Author:** Nathan Roche  
**Date:** April 2026
