# DDE Computer-Assisted Proof Project (CAPD)

This project implements computer-assisted proofs of periodic orbits for the Cubic Ikeda (a delay differential equation) using CAPD.

## Structure

- `generateImages.cpp` – main executable (visualization)
- `periodicOrbit.cpp` – main executable (proof computation)
- `computations/` – shared implementation (.cpp + .h together)
  - `pseudospectral.h/.cpp` – Chebyshev pseudospectral approximation of the delay term: Chebyshev nodes and the approximation matrix turning the DDE into a finite-dimensional ODE.
  - `attractor.h/.cpp` – Computation of the approximation of the attractor in the system and gnuplot-based visualisation of the first and last coordinates.
  - `newton.h/.cpp` – Newton's method for fixed points of the Poincaré map (`getZero`), and `getCandidate`, which continues a fixed point along a parameter path from a less chaotic to a more chaotic value of `a`.
  - `pca.h/.cpp` – Principal component analysis and plotting utilities for the attractor in new coordinates: `mainDirections`, `samplePoints`/`iteratePoincare`, and gnuplot-based 3D visualization (`plotAttractor`).
  - `cover.h/.cpp` – Covering relation checks for the computer-assisted proof (rigorous interval arithmetic verification of set coverings under the Poincaré map).
  - `bifurcation.h/.cpp` – Bifurcation diagram computation for the Cubic Ikeda map over a range of parameter values.
- `output/` – Data generated during program execution 

## Dependencies

- C++17 compiler
- CAPD library (included as submodule)
- Eigen (installed automatically by script or system package)

## Quick Start

To download and run the code it is enough to do the following in bash on Linux:
```bash
git clone --recurse-submodules https://github.com/dar-ser/Cubic-Ikeda.git
cd Cubic-Ikeda
bash setup.sh
make run
```

## More Options
For generating images of the approximated attractor from `attractor.cpp` please uncomment generateImage in `makefile` (line 2 and line 42).
For also generating the bifurcation diagram please uncomment lines 99-100 in `generateImages.cpp`.
