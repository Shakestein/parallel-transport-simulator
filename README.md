# Parallel Transport Simulator

A zero-dependency, header-only C++17 differential geometry engine that simulates and visualizes the parallel transport of tangent vectors along curves on curved 2D Riemannian manifolds embedded in $\mathbb{R}^3$.

![C++17](https://img.shields.io/badge/language-C%2B%2B17-blue.svg)
![Build](https://img.shields.io/badge/build-CMake%20%2F%20g%2B%2B-green.svg)
![License](https://img.shields.io/badge/license-MIT-informational.svg)

---

## Overview

Parallel transport is a fundamental concept in differential geometry and general relativity. On a curved surface, moving a vector along a path while keeping it "locally parallel" causes the vector to rotate relative to its starting orientation when traversing a closed loop—a geometric phenomenon known as **holonomy**.

This project provides:
- **Exact Analytical Geometry:** Surfaces (Sphere, Torus, Hyperbolic Paraboloid, Monkey Saddle) with analytically derived Christoffel symbols $\Gamma^i_{jk}$.
- **High-Precision Numerical Solver:** 4th-Order Runge-Kutta (RK4) integration of the geodesic/parallel transport differential equations with Riemannian metric norm renormalization (achieving machine-epsilon precision $\sim 10^{-13}$ rad).
- **Interactive 3D Visualizations:** Generates standalone, interactive HTML pages powered by [Plotly.js](https://plotly.com/javascript/) showing the manifold surface, the trajectory, and color-coded vector fields along the transport path.
- **Zero External Dependencies:** Built entirely with standard modern C++ (C++17 standard library only).

---

## Project Structure

```text
Parallel Transport Simulator/
├── CMakeLists.txt         # CMake build configuration
├── DOCUMENTATION.md       # Exhaustive line-by-line codebase & mathematical documentation
├── include/               # Header-only core library
│   ├── vec3.h             # 3D vector arithmetic & Euclidean geometry
│   ├── surfaces.h         # Manifold abstractions & analytical Christoffel symbols
│   ├── curves.h           # Parametric curves and velocity fields
│   ├── transport.h        # RK4 ODE integrator with metric preservation
│   └── visualization.h    # Self-contained HTML + Plotly.js generator
├── src/
│   └── main.cpp           # Interactive CLI & preset validation suite
└── output/                # Generated 3D interactive HTML plots
```

---

## Quick Start

### Build with g++ (MinGW / Clang / GCC)
```bash
g++ -std=c++17 -O2 -I include src/main.cpp -o parallel_transport.exe
./parallel_transport.exe
```

### Build with CMake
```bash
cmake -B build
cmake --build build --config Release
./build/parallel_transport.exe
```

---

## Features & Presets

The simulator includes interactive presets demonstrating canonical differential geometric phenomena:
1. **Sphere Latitude Circle:** Demonstrates non-trivial holonomy matching the solid angle $\Delta \theta = 2\pi(1 - \cos\theta_0)$.
2. **Sphere Great Circle:** Confirms geodesic parallel transport with vanishing holonomy.
3. **Torus Poloidal & Toroidal Loops:** Compares transport on positive vs. negative Gaussian curvature regions.
4. **Torus $(p, q)$ Knot:** Trajectory wrapping multiple times across toroidal homology cycles.
5. **Hyperbolic Paraboloid (Saddle):** Parallel transport on a non-compact surface of negative Gaussian curvature.
6. **Monkey Saddle:** Demonstrates behavior around higher-order degenerate critical points ($z = u^3 - 3uv^2$).
7. **Custom Simulation:** Allows user-defined surfaces, initial tangent vectors, and curve configurations.

---

## Documentation

For an exhaustive, line-by-line mathematical and programming breakdown, see [DOCUMENTATION.md](DOCUMENTATION.md).

---

## License

MIT License. Feel free to use, modify, and distribute.
