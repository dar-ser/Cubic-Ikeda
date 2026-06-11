# DDE Computer-Assisted Proof Project (CAPD)

This project implements computer-assisted proofs of periodic orbits for the Cubic Ikeda (a delay differential equation) using CAPD.

## Structure

- `generateImages.cpp` – main executable (visualization)
- `periodicOrbit.cpp` – main executable (proof computation)
- `computations/` – shared implementation (.cpp + .h together, CAPD-safe)
- `external/capd/` – CAPD library (git submodule)
- `output/` – generated results

## Dependencies

- C++17 compiler
- CAPD library (included as submodule)
- Eigen (installed automatically by script or system package)

## Quick Start

```bash
git clone --recurse-submodules <repo>
cd <repo>
bash setup.sh
make run
