# Heterogeneous Phase-Field Fracture Simulation 🚀

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Language: C++](https://img.shields.io/badge/Language-C%2B%2B-blue.svg)
![Framework: deal.II](https://img.shields.io/badge/Framework-deal.II-green.svg)
![Parallel: MPI](https://img.shields.io/badge/Parallel-MPI-red.svg)

This repository contains a high-performance Phase-Field Fracture Mechanics (PFM) simulation designed to model crack propagation in heterogeneous materials. The current configuration simulates an **Aluminum-infiltrated Steel Matrix**, capturing complex crack-inclusion interactions.

## 🌟 Key Features

*   **Heterogeneous Material Support**: Distinct physical properties ($\lambda, \mu, G_c$) for Steel (Matrix) and Aluminum (Inclusions) mapped directly from Gmsh Physical Tags.
*   **Adaptive Mesh Refinement (AMR)**: Automatically refines the mesh near the crack tip and damage zones to ensure high resolution where it matters most.
*   **Robust Staggered Solver**: Implements a staggered approach (displacement-damage coupling) with a **Dynamic Tolerance Fallback** mechanism to ensure convergence during complex propagation steps.
*   **Automated Post-Processing**: Includes Python scripts for professional results visualization and report generation.

## 🖼️ Visual Results

### Force-Displacement Comparison
Comparison of the simulated results against experimental/reference data from literature.
![Force-Displacement Comparison](images/force_displacement_comparison_hetero.png)

### Material Configuration & Crack Propagation
| Initial Geometry & Inclusions | Crack Tip & Notch | Final Crack Pattern |
| :---: | :---: | :---: |
| ![Inclusions](images/inclusion%20in%20material.png) | ![Initial Crack](images/Initial_crack.png) | ![Cracked Plate](images/cracked_plate.png) |

## 📂 Project Structure

| Directory/File | Description |
| :--- | :--- |
| `source/` | C++ source files containing solver logic, assembly, and physics. |
| `include/` | Header files for the `PhaseField` class and variable definitions. |
| `mesh/` | Gmsh `.geo` scripts and generated `.msh` files. |
| `images/` | **(New)** Visual results, plots, and crack propagation snapshots. |
| `output/` | Simulation results (CSV data). |
| `report/` | LaTeX source for automated report generation. |
| `Phase_Field.cc` | The main entry point and runtime configuration. |
| `plot_force_displacement_comparison.py` | Professional Matplotlib script for results visualization. |

## 🛠️ Installation & Prerequisites

### Prerequisites
*   **deal.II** (9.x recommended)
*   **Gmsh** (for mesh generation)
*   **MPI** (for parallel execution)
*   **Python 3.x** (with `matplotlib`)

### Build Instructions
```bash
# Generate the mesh
gmsh -2 mesh/geometry.geo -o mesh/geometry.msh

# Compile the project
cmake .
make -j$(nproc)
```

## 🚀 Usage

Run the simulation:
```bash
./excecutable_output
```

Generate the comparison plot manually:
```bash
python3 plot_force_displacement_comparison.py --generated output/force_displacement.csv --validation force_displacement_VALIDATION.ods
```

## 📈 Numerical Formulation

The simulation solves the coupled system for displacement $\mathbf{u}$ and damage field $\phi$:
1.  **Elastic energy splitting**: Uses spectral decomposition to prevent crack growth under pure compression.
2.  **Staggered Scheme**: 
    *   Solve for $\mathbf{u}$ (Fixed $\phi$).
    *   Solve for $\phi$ (Fixed $\mathbf{u}$).
    *   Iterate until $L_2$ norm of the update is below tolerance.

---
**Author**: [KALPESH SINGH (kalp-dmech)](https://github.com/kalp-dmech) | [LinkedIn](https://www.linkedin.com/in/kalpeshiitm/)
**Research Area**: Computational Fracture Mechanics, Phase-Field Methods.
