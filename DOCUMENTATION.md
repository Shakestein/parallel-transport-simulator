# Parallel Transport Simulator — Complete Documentation

This document provides a comprehensive guide to the Parallel Transport Simulator, covering the project's scaffolding, a line-by-line explanation of the mathematical and programming logic, and a detailed user guide.

---

## 1. Project Scaffolding and Architecture

```text
Parallel Transport Simulator/
├── CMakeLists.txt                  ← Build system configuration
├── parallel_transport.exe          ← Compiled binary
├── include/                        ← Header files (core logic)
│   ├── vec3.h                      ← 3D vector math primitives
│   ├── surfaces.h                  ← Surface definitions (Sphere, Torus, etc.)
│   ├── curves.h                    ← Parametric curve definitions
│   ├── transport.h                 ← RK4 numerical ODE solver
│   └── visualization.h             ← HTML + Plotly.js code generator
├── src/
│   └── main.cpp                    ← Application entry point & CLI
└── output/                         ← Generated HTML visualizations
```

### Why a Header-Only Library?
The core simulation engine (`include/`) is designed as a **header-only library**. This means there are no `.cpp` files for the core classes. 
- **Simplicity:** You only need to compile `main.cpp`. The compiler pulls in everything else automatically.
- **Zero Dependencies:** The project relies solely on the C++17 Standard Library. No external math or graphics libraries (like Eigen or OpenGL) are required.
- **Inlining:** Compilers can aggressively inline functions defined in headers, which is crucial for the performance of the RK4 integrator running thousands of steps.

---

## 2. Line-by-Line Code Explanation

### 2.1 `vec3.h` (Vector Mathematics)
This file defines a lightweight 3D vector struct used for spatial calculations.

- **Lines 1-4:** `#pragma once` prevents double-inclusion. We include `<cmath>` for square roots and `<iostream>` for console output.
- **Lines 6-10:** `struct Vec3` defines the vector. Members `x, y, z` are `double` for high precision (crucial for numerical integration). We provide a default constructor (initializes to 0) and a parameterized constructor.
- **Lines 12-19:** Operator overloading (`+`, `-`, `*`, `/` and their `+=` variants). This allows natural mathematical syntax like `v1 + v2 * 0.5`. All scalar operations are per-component.
- **Line 21:** `dot(const Vec3& o)` computes the Euclidean dot product ($x_1x_2 + y_1y_2 + z_1z_2$). Essential for calculating lengths, angles, and the surface metric.
- **Lines 23-27:** `cross(const Vec3& o)` computes the cross product. Used to calculate surface normals ($\mathbf{r}_u \times \mathbf{r}_v$).
- **Lines 29-36:** `norm()` and `norm_sq()` calculate the vector's magnitude. `normalized()` returns a unit vector, including a safety check `if (n < 1e-15)` to prevent division by zero.
- **Lines 38-43:** `friend` functions. One allows scalar multiplication from the left (`2.0 * v`), and the other overloads `<<` to easily print vectors to `std::cout`.

### 2.2 `surfaces.h` (Differential Geometry)
Defines the mathematical manifolds.

- **Lines 17-77:** `class Surface`. This is an abstract base class (interface).
  - `virtual ~Surface() = default;` ensures proper cleanup of derived classes.
  - Pure virtual functions (`= 0`) force derived classes to define their `name`, `position` $\mathbf{r}(u,v)$, partial derivatives `r_u` and `r_v`, `christoffel` symbols ($\Gamma^i_{jk}$), and parameter `domain`.
  - `normal()` computes the unit normal $\hat{\mathbf{n}} = \frac{\mathbf{r}_u \times \mathbf{r}_v}{|\mathbf{r}_u \times \mathbf{r}_v|}$.
  - `metric()` computes the First Fundamental Form coefficients: $E = \mathbf{r}_u \cdot \mathbf{r}_u$, $F = \mathbf{r}_u \cdot \mathbf{r}_v$, $G = \mathbf{r}_v \cdot \mathbf{r}_v$.
  - `metric_norm()` computes the length of an intrinsic tangent vector $(V^1, V^2)$ using the metric tensor: $\sqrt{E(V^1)^2 + 2FV^1V^2 + G(V^2)^2}$.
- **Lines 79-132:** `class Sphere`. Parametrized by colatitude $u \in (0, \pi)$ and longitude $v \in (0, 2\pi)$. 
  - The domain excludes exact poles (`0.05` padding) to avoid metric singularities where $G=0$.
  - `christoffel()` implements the non-zero symbols for a sphere: $\Gamma^\theta_{\phi\phi} = -\sin\theta\cos\theta$ and $\Gamma^\phi_{\theta\phi} = \cot\theta$.
- **Lines 134-186:** `class Torus`. Parametrized by poloidal angle $u$ and toroidal angle $v$. Christoffel symbols are derived analytically.
- **Lines 188-226:** `class HyperbolicParaboloid`. A saddle surface defined as a Monge patch $z = (u^2 - v^2)/2$. The symbols are calculated using the general Monge patch formula $\Gamma^i_{jk} = \frac{h_i h_{jk}}{1 + h_u^2 + h_v^2}$.
- **Lines 228-279:** `class MonkeySaddle`. $z = u^3 - 3uv^2$. Calculates symbols directly using analytical partial derivatives of the height function $h(u,v)$.

### 2.3 `curves.h` (Parametric Paths)
Defines paths $(u(t), v(t))$ on the surfaces.

- **Lines 12-15:** `clamp_val()` is a custom fallback for `std::clamp` to support older MinGW compilers.
- **Lines 20-33:** `class Curve` abstract base class. `evaluate()` must provide position $(u,v)$ and velocity $(\dot{u}, \dot{v})$ for a given parameter $t$.
- **Lines 35-59:** `LatitudeCircle`. Traces a horizontal circle on a sphere. $u(t) = \theta_0$ (constant), $v(t) = t$. Velocity is $(0, 1)$.
- **Lines 61-108:** `GreatCircle`. A geodesic on a sphere. Because converting a tilted Cartesian circle to spherical coordinates is complex, this class uses central finite differences (`dt = 1e-7`) to numerically calculate $\dot{u}$ and $\dot{v}$. It includes wrap-around handling for the azimuth angle $\phi$.
- **Lines 110-211:** Other curves (`SpiralCurve`, `CircularPath`, `TorusKnot`, `ToroidalLoop`) implement specific $u(t)$ and $v(t)$ equations and their exact analytical derivatives.

### 2.4 `transport.h` (The RK4 Solver)
The computational heart of the simulator. Integrates the parallel transport ODE.

- **Lines 23-37:** `TransportResult` struct. Stores arrays of all states at every time step (positions, vectors, norms) so they can be visualized later. Also stores error metrics.
- **Lines 42-132:** `ParallelTransportSolver::solve()`. 
  - Calculates step size `dt`.
  - Reserves memory in `std::vector`s to prevent reallocation overhead.
  - Enters the `for` loop for `n_steps`.
  - **RK4 Integration:** Evaluates the derivative `transport_rhs` at 4 points (current `k1`, two midpoints `k2, k3`, and next point `k4`). Averages them using the standard Runge-Kutta weights $(1/6, 1/3, 1/3, 1/6)$ to update the intrinsic vector $V^1, V^2$.
  - **Metric Renormalization:** Due to floating-point drift, the vector's length might slightly change. Lines 117-122 project the vector back to its exact initial `initial_norm` using scalar multiplication, preserving the Levi-Civita metric compatibility.
- **Lines 134-156:** `transport_rhs()`. Implements the ODE: $\frac{dV^i}{dt} = -\Gamma^i_{jk} V^j \frac{du^k}{dt}$. It iterates over $j,k \in \{0,1\}$ and accumulates the sum.
- **Lines 177-215:** `compute_diagnostics()`. Calculates the holonomy angle using $\arccos(\mathbf{V}_{start} \cdot \mathbf{V}_{end})$. Calculates the maximum drift from the initial norm, and the maximum tangency error (deviation from the tangent plane, which should be $\approx 0$).

### 2.5 `visualization.h` (HTML & Plotly Generation)
Generates an interactive 3D web view.

- **Lines 20-408:** `Visualization::generate()`.
  - Creates a 2D mesh grid over the surface's $(u,v)$ domain.
  - Constructs a massive `std::ostringstream html` containing a complete HTML document.
  - Embeds CSS for a sleek dark mode UI.
  - Embeds Plotly.js via CDN.
  - Dumps the C++ vectors (surface mesh, curve points, and arrow vectors) into JavaScript arrays using JSON formatting (`to_json_2d`, `to_json_1d_vec3`).
  - Creates a `surface` trace, a `scatter3d` line trace for the curve, and generates 3D arrows (a line for the shaft, a `cone` for the tip) at sampled intervals (`arrow_skip`).
  - Colors the first arrow Cyan, the last Orange, and gradients the middle arrows.
  - Uses system commands (`mkdir`) to ensure the `output/` directory exists and writes the file.
- **Lines 410-422:** `open_in_browser()`. Uses OS-specific shell commands (`start`, `open`, `xdg-open`) to automatically open the resulting HTML file.

### 2.6 `main.cpp` (Application CLI)
The entry point connecting everything.

- **Lines 20-76:** `run_simulation()`. Wrapper function that calls the solver, prints the diagnostics cleanly to the console (comparing expected vs. computed holonomy), and triggers the visualization.
- **Lines 78-150:** Preset demos. Functions like `demo_sphere_latitude()` instantiate specific surfaces and curves. It calculates the theoretical expected holonomy. For a latitude circle at colatitude $\theta$, the enclosed solid angle (and thus holonomy) is $2\pi(1 - \cos\theta)$.
- **Lines 152-281:** `custom_simulation()`. Provides an interactive wizard via `std::cin` to let the user mix and match surfaces, curves, initial vectors, and step counts.
- **Lines 283-340:** `main()`. Prints the ASCII banner and runs the infinite menu loop.

---

## 3. User Guide

### Requirements
- A modern C++ compiler (GCC, Clang, or MSVC) supporting C++17.
- An internet connection when viewing results (to load the Plotly.js library from CDN).

### How to Build

**Option A: Direct Compilation (Recommended)**
Open a terminal in the project root and run:
```bash
g++ -std=c++17 -O2 -I include src/main.cpp -o parallel_transport.exe
```

**Option B: CMake**
```bash
cmake -B build
cmake --build build --config Release
```

### How to Run
Execute the binary:
```bash
./parallel_transport.exe
```

### Using the Application
Upon launching, you will be presented with a menu:
1. **Presets (Options 1-7):** Choose a preset to instantly run a predefined simulation. 
   - *Example:* Option 1 runs a vector around a 45-degree latitude circle on a sphere. This demonstrates non-zero holonomy (the vector points in a different direction when it returns).
   - *Example:* Option 3 runs a vector along a Great Circle. Because great circles are geodesics, the holonomy is exactly 0.
2. **Custom Simulation (Option 8):** Build your own setup.
   - Choose a surface (e.g., Torus).
   - Choose parameters for that surface (e.g., Major/Minor radius).
   - Choose a curve (e.g., Toroidal Loop).
   - Specify the initial vector components in the intrinsic $u, v$ basis. `(1, 0)` means pointing exactly along the $u$ parameter direction.
   - Set the number of integration steps (default 2000 is highly accurate).
3. **Viewing Results:**
   - The CLI will print the math results, showing the Holonomy Angle, Norm Drift (error), and Tangency Error.
   - Your default web browser will automatically open showing an interactive 3D plot.
   - You can click and drag to rotate the 3D plot, scroll to zoom, and hover over elements to see exact coordinate data.

### Troubleshooting
- **Browser doesn't open?** Look inside the `output/` directory in the project folder and double-click the newly generated `.html` file.
- **Blank web page?** Ensure you have an active internet connection so the browser can download the Plotly graphics library.
