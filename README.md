# hyphal-flow

A computational framework for studying a drop flowing through a single fungal hypha branch with viscoelastic effects. Built on Basilisk C, this code simulates three-phase non-coalescing systems, allowing you to track fluid-structure interactions among the drop, hyphal wall (treated as a viscoelastic Kelvin–Voigt solid), and surrounding medium (cytoplasm).

## Key Features

### Three-Phase Model with Viscoelasticity
- Drop–Hypha–Cytoplasm phases with separate density, viscosity, and elastic properties
- Log-conformation viscoelastic model to handle elasticity in the hypha or fluid
- Surface tension and non-coalescing mechanics between phases

### Adaptive Mesh Refinement
- Automated refinement based on:
  - Phase fraction field
  - Interfacial curvature
  - Velocity and stress gradients

### High-Performance Capability
- Compatible with OpenMP for shared-memory parallelism
- Optional MPI support for distributed-memory HPC
- Customizable event-driven time stepping for transient simulations

### Configurable Physics
- User-defined non-dimensional groups (Ohnesorge, Bond, Deborah, Elasto-capillary)
- Flexible boundary conditions and domain setup
- Explicit or implicit viscous terms based on Basilisk's Navier–Stokes solvers

## Project Structure

```
hyphal-flow/
├── .vscode/                              # VS Code configuration
├── basilisk/                             # Core Basilisk C source and headers
├── postProcess/
│   ├── getData-elastic-nonCoalescence.c  # Data extraction for visualization
│   └── getFacet-threePhase.c            # Interface extraction utility
├── src-local/
│   ├── log-conform-elastic.h             # Log-conformation for elastic fluids
│   ├── log-conform-viscoelastic.h        # Extended viscoelastic model
│   ├── reduced-three-phase-nonCoalescing.h
│   ├── three-phase-nonCoalescing-elastic.h
│   └── three-phase-nonCoalescing-viscoelastic.h
├── testCases/
│   ├── hypha.c                           # Main simulation file
│   ├── Makefile                          # Compilation instructions
│   └── runCodesInParallel.sh            # MPI parallel run script
├── LICENSE
├── README.md
└── reset_install_requirements.sh         # Basilisk environment setup
```

## Installation

### Prerequisites
- A C compiler (e.g., GCC 7.0+)
- OpenMP (often bundled with GCC)
- MPI implementation (e.g., OpenMPI or MPICH) for parallel execution
- Basilisk C source code (managed through this repository)

### Setup Steps
1. Clone the repository:
```bash
git clone https://github.com/VatsalSy/hyphal-flow.git
cd hyphal-flow
```

2. Set up the environment:
```bash
./reset_install_requirements.sh
```
This script:
- Installs or updates Basilisk in a local directory
- Configures environment variables in .project_config

3. Verify your setup:
```bash
source .project_config
qcc --version
```

## Usage

### Compiling the Code

#### Single-Core / OpenMP
```bash
qcc -Wall -O2 -fopenmp -I$(PWD)/src-local testCases/hypha.c -o hypha -lm

# Run with desired threads
export OMP_NUM_THREADS=4
./hypha
```

#### MPI (Distributed Memory)
```bash
CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 \
    -I$(PWD)/src-local testCases/hypha.c -o hyphaMPI -lm
mpirun -np 8 ./hyphaMPI
```

### Simulation Parameters

In `hypha.c`, you can adjust:
- Mesh resolution (MAXlevel / MINlevel)
- Non-dimensional numbers (Ohd, Ohf, Ohc, Ec_*, De_*)
- Bond number and domain size (Ldomain)
- Time stepping and output frequency

### Post-Processing

The `postProcess` folder contains utilities for data analysis:
- `getData-elastic-nonCoalescence.c`: Extracts fields for visualization
- `getFacet-threePhase.c`: Retrieves interface geometry

## Contributing

Contributions are welcome! To propose changes:
1. Fork the repository
2. Create a feature branch: `git checkout -b feature/MyFeature`
3. Commit and push your changes
4. Open a Pull Request

Please ensure new features are tested and documented.

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this code in your research, please cite:

```bibtex
@software{hyphal_flow_2024,
  author       = {Vatsal Sanjay},
  title        = {Hyphal Flow: A Three-Phase Viscoelastic Framework},
  year         = {2024},
  version      = {v1.0},
  url          = {https://github.com/VatsalSy/hyphal-flow}
}
```

## Acknowledgments
- Built upon the Basilisk C framework by Stéphane Popinet.
- Special thanks to the Physics of Fluids group for support. 
- Special thanks to Mazi Jalaal and Eric Lauga for insights on fungal fluid dynamics.
