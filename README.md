# Computational Physics Simulations in C++

This repository contains a collection of common numerical simulations and algorithms implemented in C++. It is actively maintained as part of a university course, and more projects will be added to this repository as the semester progresses.

## Repository Structure

The repository is organized as a collection of independent sub-projects. Currently, it includes the following simulation directories:

* **`Jacobi_Rotation_Algorithm/`**: Numerical methods for solving eigenvalue problems[cite: 10].
* **`Penning_Trap/`**: Object-oriented simulation of particle dynamics within electromagnetic fields[cite: 10].

### Per-Project Configuration

Because each simulation addresses a different physical system and numerical method, every folder is entirely self-contained. Inside each individual project directory, you will find:
* A dedicated **`README.md`** detailing the specific physical theory, mathematical derivations, and usage instructions.
* Project-specific **requirements** and build configurations (e.g., Makefiles or CMake configurations).
* A local **`.clang-format`** file to enforce consistent C++ code styling. 

## Getting Started

To explore or run a specific simulation, navigate to its respective folder and follow the build and execution instructions provided in its local `README.md`.

```bash
# Example workflow
git clone <your-repository-url>
cd Penning_Trap
# Follow the instructions in Penning_Trap/README.md
