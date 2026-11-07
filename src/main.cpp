#include "vec3.h"
#include "surfaces.h"
#include "curves.h"
#include "transport.h"
#include "visualization.h"

#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <iomanip>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ============================================================================
// Run a simulation and display results
// ============================================================================
void run_simulation(const Surface& surface, const Curve& curve,
                    double V1_0, double V2_0, int n_steps,
                    const std::string& description = "",
                    double expected_holonomy = -1.0) {
    std::cout << "\n";
    std::cout << "========================================================\n";
    if (!description.empty()) {
        std::cout << "  " << description << "\n";
    }
    std::cout << "  Surface: " << surface.name() << "\n";
    std::cout << "  Curve:   " << curve.name() << "\n";
    std::cout << "  Initial vector: (" << V1_0 << ", " << V2_0 << ")\n";
    std::cout << "  Steps:   " << n_steps << "\n";
    std::cout << "========================================================\n";

    std::cout << "  Integrating parallel transport equation (RK4)...\n";
    TransportResult result = ParallelTransportSolver::solve(
        surface, curve, V1_0, V2_0, n_steps
    );

    // Print diagnostics
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\n  --- Results ---\n";
    std::cout << "  Holonomy angle:     " << result.holonomy_angle << " rad ("
              << result.holonomy_angle * 180.0 / M_PI << " deg)\n";

    if (expected_holonomy >= 0.0) {
        double error = std::abs(result.holonomy_angle - expected_holonomy);
        std::cout << "  Expected holonomy:  " << expected_holonomy << " rad ("
                  << expected_holonomy * 180.0 / M_PI << " deg)\n";
        std::cout << "  Holonomy error:     " << std::scientific << error << std::fixed << "\n";
    }

    std::cout << std::scientific << std::setprecision(3);
    std::cout << "  Max norm drift:     " << result.max_norm_drift
              << " (should be ~0)\n";
    std::cout << "  Max tangency error: " << result.max_tangency_error
              << " (should be ~0)\n";
    std::cout << std::fixed << std::setprecision(6);

    // Initial and final 3D vectors
    Vec3 v_start = result.vectors_3d.front();
    Vec3 v_end   = result.vectors_3d.back();
    std::cout << "  Initial 3D vector:  " << v_start << "  |V| = " << v_start.norm() << "\n";
    std::cout << "  Final 3D vector:    " << v_end   << "  |V| = " << v_end.norm()   << "\n";

    // Generate visualization
    std::cout << "\n  Generating visualization...\n";
    std::string filepath = Visualization::generate(surface, curve, result);

    if (!filepath.empty()) {
        std::cout << "  Opening in browser...\n";
        Visualization::open_in_browser(filepath);
    }

    std::cout << "========================================================\n\n";
}

// ============================================================================
// Preset demos
// ============================================================================
void demo_sphere_latitude() {
    Sphere sphere(1.0);
    double theta0 = M_PI / 4.0;  // 45° colatitude
    LatitudeCircle curve(theta0);

    // Initial vector: unit vector along ∂/∂θ direction
    // V = (1, 0) means V = ∂r/∂θ
    // The holonomy (rotation angle) for a latitude circle at colatitude θ
    // is the solid angle of the enclosed spherical cap = 2π(1 - cosθ)
    // This equals what acos(V_start · V_end) measures
    double expected = 2.0 * M_PI * (1.0 - cos(theta0));

    run_simulation(sphere, curve, 1.0, 0.0, 2000,
                   "HOLONOMY DEMO: Vector transport around latitude circle on sphere",
                   expected);
}

void demo_sphere_latitude_60() {
    Sphere sphere(1.0);
    double theta0 = M_PI / 3.0;  // 60° colatitude
    LatitudeCircle curve(theta0);

    double expected = 2.0 * M_PI * (1.0 - cos(theta0));

    run_simulation(sphere, curve, 1.0, 0.0, 2000,
                   "Sphere: latitude circle at 60 deg colatitude",
                   expected);
}

void demo_sphere_great_circle() {
    Sphere sphere(1.0);
    GreatCircle curve(M_PI / 4.0);  // tilted 45°

    // On a geodesic, parallel transport of the tangent vector results in zero holonomy
    run_simulation(sphere, curve, 1.0, 0.0, 2000,
                   "Sphere: great circle (geodesic) — expect zero holonomy",
                   0.0);
}

void demo_torus_toroidal() {
    Torus torus(2.0, 0.7);
    ToroidalLoop curve(0.0);  // outer equator (u=0)

    run_simulation(torus, curve, 1.0, 0.0, 2000,
                   "Torus: toroidal loop along outer equator");
}

void demo_torus_spiral() {
    Torus torus(2.0, 0.7);
    SpiralCurve curve(0.0, 0.5, 1.0, 2.0 * M_PI);  // gentle spiral

    run_simulation(torus, curve, 1.0, 0.0, 2000,
                   "Torus: spiral path");
}

void demo_saddle_circle() {
    HyperbolicParaboloid saddle;
    CircularPath curve(0.0, 0.0, 0.8);

    run_simulation(saddle, curve, 1.0, 0.0, 2000,
                   "Hyperbolic Paraboloid: circular path around origin");
}

void demo_monkey_saddle_circle() {
    MonkeySaddle monkey;
    CircularPath curve(0.0, 0.0, 0.7);

    run_simulation(monkey, curve, 1.0, 0.0, 2000,
                   "Monkey Saddle: circular path around origin");
}

void custom_simulation() {
    std::cout << "\n=== Custom Simulation Setup ===\n\n";

    // Choose surface
    std::cout << "Choose surface:\n";
    std::cout << "  1. Sphere\n";
    std::cout << "  2. Torus\n";
    std::cout << "  3. Hyperbolic Paraboloid\n";
    std::cout << "  4. Monkey Saddle\n";
    std::cout << "  > ";
    int surf_choice;
    std::cin >> surf_choice;

    std::unique_ptr<Surface> surface;
    switch (surf_choice) {
        case 1: {
            double R;
            std::cout << "  Sphere radius [1.0]: ";
            std::string input;
            std::getline(std::cin >> std::ws, input);
            R = input.empty() ? 1.0 : std::stod(input);
            surface = std::make_unique<Sphere>(R);
            break;
        }
        case 2: {
            double R, r;
            std::cout << "  Major radius [2.0]: ";
            std::string input;
            std::getline(std::cin >> std::ws, input);
            R = input.empty() ? 2.0 : std::stod(input);
            std::cout << "  Minor radius [0.7]: ";
            std::getline(std::cin, input);
            r = input.empty() ? 0.7 : std::stod(input);
            surface = std::make_unique<Torus>(R, r);
            break;
        }
        case 3:
            surface = std::make_unique<HyperbolicParaboloid>();
            break;
        case 4:
            surface = std::make_unique<MonkeySaddle>();
            break;
        default:
            std::cout << "Invalid choice.\n";
            return;
    }

    // Choose curve
    std::cout << "\nChoose curve:\n";
    std::cout << "  1. Latitude Circle (for sphere)\n";
    std::cout << "  2. Great Circle (for sphere)\n";
    std::cout << "  3. Toroidal Loop (for torus)\n";
    std::cout << "  4. Spiral\n";
    std::cout << "  5. Circular Path in parameter space\n";
    std::cout << "  > ";
    int curve_choice;
    std::cin >> curve_choice;

    std::unique_ptr<Curve> curve;
    switch (curve_choice) {
        case 1: {
            double theta;
            std::cout << "  Colatitude (radians) [0.785]: ";
            std::string input;
            std::getline(std::cin >> std::ws, input);
            theta = input.empty() ? M_PI / 4.0 : std::stod(input);
            curve = std::make_unique<LatitudeCircle>(theta);
            break;
        }
        case 2: {
            double tilt;
            std::cout << "  Tilt angle (radians) [0.785]: ";
            std::string input;
            std::getline(std::cin >> std::ws, input);
            tilt = input.empty() ? M_PI / 4.0 : std::stod(input);
            curve = std::make_unique<GreatCircle>(tilt);
            break;
        }
        case 3: {
            double u0;
            std::cout << "  Poloidal angle [0.0]: ";
            std::string input;
            std::getline(std::cin >> std::ws, input);
            u0 = input.empty() ? 0.0 : std::stod(input);
            curve = std::make_unique<ToroidalLoop>(u0);
            break;
        }
        case 4: {
            curve = std::make_unique<SpiralCurve>();
            break;
        }
        case 5: {
            double cu, cv, r;
            std::cout << "  Center u [0.0]: ";
            std::string input;
            std::getline(std::cin >> std::ws, input);
            cu = input.empty() ? 0.0 : std::stod(input);
            std::cout << "  Center v [0.0]: ";
            std::getline(std::cin, input);
            cv = input.empty() ? 0.0 : std::stod(input);
            std::cout << "  Radius [0.8]: ";
            std::getline(std::cin, input);
            r = input.empty() ? 0.8 : std::stod(input);
            curve = std::make_unique<CircularPath>(cu, cv, r);
            break;
        }
        default:
            std::cout << "Invalid choice.\n";
            return;
    }

    // Initial vector
    double V1, V2;
    std::cout << "\nInitial vector components (V1, V2):\n";
    std::cout << "  V1 [1.0]: ";
    std::string input;
    std::getline(std::cin >> std::ws, input);
    V1 = input.empty() ? 1.0 : std::stod(input);
    std::cout << "  V2 [0.0]: ";
    std::getline(std::cin, input);
    V2 = input.empty() ? 0.0 : std::stod(input);

    // Number of steps
    int n_steps;
    std::cout << "\nNumber of integration steps [2000]: ";
    std::getline(std::cin, input);
    n_steps = input.empty() ? 2000 : std::stoi(input);

    run_simulation(*surface, *curve, V1, V2, n_steps, "Custom Simulation");
}

// ============================================================================
// Main
// ============================================================================
int main() {
    std::cout << R"(
  ____                 _ _      _   _____                                _
 |  _ \ __ _ _ __ __ _| | | ___| | |_   _| __ __ _ _ __  ___ _ __   ___| |_
 | |_) / _` | '__/ _` | | |/ _ \ |   | || '__/ _` | '_ \/ __| '_ \ / _ \ __|
 |  __/ (_| | | | (_| | | |  __/ |   | || | | (_| | | | \__ \ |_) | (_) | |_
 |_|   \__,_|_|  \__,_|_|_|\___|_|   |_||_|  \__,_|_| |_|___/ .__/ \___/ \__|
                                                              |_|
               S I M U L A T O R   on   C U R V E D   S U R F A C E S
)" << "\n";

    while (true) {
        std::cout << "Choose a demo or build a custom simulation:\n\n";
        std::cout << "  --- Sphere ---\n";
        std::cout << "  1. Latitude circle at 45 deg (classic holonomy demo)\n";
        std::cout << "  2. Latitude circle at 60 deg\n";
        std::cout << "  3. Great circle (geodesic — zero holonomy)\n\n";
        std::cout << "  --- Torus ---\n";
        std::cout << "  4. Toroidal loop (outer equator)\n";
        std::cout << "  5. Spiral path\n\n";
        std::cout << "  --- Saddle Surfaces ---\n";
        std::cout << "  6. Hyperbolic paraboloid — circular path\n";
        std::cout << "  7. Monkey saddle — circular path\n\n";
        std::cout << "  --- Custom ---\n";
        std::cout << "  8. Custom simulation (choose surface, curve, vector)\n\n";
        std::cout << "  0. Exit\n\n";
        std::cout << "  > ";

        int choice;
        if (!(std::cin >> choice)) {
            std::cin.clear();
            std::cin.ignore(10000, '\n');
            continue;
        }

        switch (choice) {
            case 0:
                std::cout << "\nGoodbye!\n";
                return 0;
            case 1: demo_sphere_latitude();       break;
            case 2: demo_sphere_latitude_60();    break;
            case 3: demo_sphere_great_circle();   break;
            case 4: demo_torus_toroidal();        break;
            case 5: demo_torus_spiral();          break;
            case 6: demo_saddle_circle();         break;
            case 7: demo_monkey_saddle_circle();  break;
            case 8: custom_simulation();          break;
            default:
                std::cout << "Invalid choice. Try again.\n\n";
                break;
        }
    }

    return 0;
}
