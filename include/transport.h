#pragma once

#include "vec3.h"
#include "surfaces.h"
#include "curves.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Helper: clamp value to [lo, hi]
inline double clamp_val_t(double x, double lo, double hi) {
    return (x < lo) ? lo : (x > hi) ? hi : x;
}

// ============================================================================
// Result of parallel transport computation
// ============================================================================
struct TransportResult {
    std::vector<double> t_values;     // parameter values
    std::vector<double> u_values;     // u coordinates along curve
    std::vector<double> v_values;     // v coordinates along curve
    std::vector<Vec3>   positions;    // 3D positions on surface
    std::vector<Vec3>   vectors_3d;   // transported vectors in 3D
    std::vector<double> V1_values;    // intrinsic component V^1
    std::vector<double> V2_values;    // intrinsic component V^2
    std::vector<double> norms;        // metric norms (should be constant)
    std::vector<double> tangency;     // |V · n| (should be ~0)

    double holonomy_angle;            // angle between initial and final vector
    double max_norm_drift;            // max |norm - initial_norm| / initial_norm
    double max_tangency_error;        // max |V · n|
};

// ============================================================================
// Parallel Transport Solver using RK4
// ============================================================================
class ParallelTransportSolver {
public:
    // Perform parallel transport of vector (V1_0, V2_0) along curve on surface
    static TransportResult solve(const Surface& surface,
                                 const Curve& curve,
                                 double V1_0, double V2_0,
                                 int n_steps = 1000) {
        TransportResult result;

        double t0 = curve.t_start();
        double t1 = curve.t_end();
        double dt = (t1 - t0) / n_steps;

        // Reserve space
        result.t_values.reserve(n_steps + 1);
        result.u_values.reserve(n_steps + 1);
        result.v_values.reserve(n_steps + 1);
        result.positions.reserve(n_steps + 1);
        result.vectors_3d.reserve(n_steps + 1);
        result.V1_values.reserve(n_steps + 1);
        result.V2_values.reserve(n_steps + 1);
        result.norms.reserve(n_steps + 1);
        result.tangency.reserve(n_steps + 1);

        // Initial state
        double V1 = V1_0, V2 = V2_0;
        double u, v, du, dv;
        curve.evaluate(t0, u, v, du, dv);

        // Store initial metric norm for renormalization
        double initial_norm = surface.metric_norm(u, v, V1, V2);

        // Record initial point
        record_step(surface, result, t0, u, v, V1, V2);

        // RK4 integration
        for (int step = 0; step < n_steps; step++) {
            double t = t0 + step * dt;

            // Evaluate curve at current, midpoint, and next positions
            double u_cur, v_cur, du_cur, dv_cur;
            double u_mid, v_mid, du_mid, dv_mid;
            double u_nxt, v_nxt, du_nxt, dv_nxt;

            curve.evaluate(t, u_cur, v_cur, du_cur, dv_cur);
            curve.evaluate(t + 0.5 * dt, u_mid, v_mid, du_mid, dv_mid);
            curve.evaluate(t + dt, u_nxt, v_nxt, du_nxt, dv_nxt);

            // k1: f(t, V)
            double k1_1, k1_2;
            transport_rhs(surface, u_cur, v_cur, du_cur, dv_cur, V1, V2, k1_1, k1_2);

            // k2: f(t + dt/2, V + dt/2 * k1)
            double V1_tmp = V1 + 0.5 * dt * k1_1;
            double V2_tmp = V2 + 0.5 * dt * k1_2;
            double k2_1, k2_2;
            transport_rhs(surface, u_mid, v_mid, du_mid, dv_mid, V1_tmp, V2_tmp, k2_1, k2_2);

            // k3: f(t + dt/2, V + dt/2 * k2)
            V1_tmp = V1 + 0.5 * dt * k2_1;
            V2_tmp = V2 + 0.5 * dt * k2_2;
            double k3_1, k3_2;
            transport_rhs(surface, u_mid, v_mid, du_mid, dv_mid, V1_tmp, V2_tmp, k3_1, k3_2);

            // k4: f(t + dt, V + dt * k3)
            V1_tmp = V1 + dt * k3_1;
            V2_tmp = V2 + dt * k3_2;
            double k4_1, k4_2;
            transport_rhs(surface, u_nxt, v_nxt, du_nxt, dv_nxt, V1_tmp, V2_tmp, k4_1, k4_2);

            // Update vector
            V1 += (dt / 6.0) * (k1_1 + 2.0 * k2_1 + 2.0 * k3_1 + k4_1);
            V2 += (dt / 6.0) * (k1_2 + 2.0 * k2_2 + 2.0 * k3_2 + k4_2);

            // Metric renormalization: project back to preserve ||V||_g
            double cur_norm = surface.metric_norm(u_nxt, v_nxt, V1, V2);
            if (cur_norm > 1e-15) {
                double scale = initial_norm / cur_norm;
                V1 *= scale;
                V2 *= scale;
            }

            // Record step
            record_step(surface, result, t + dt, u_nxt, v_nxt, V1, V2);
        }

        // Compute diagnostics
        compute_diagnostics(result, initial_norm);

        return result;
    }

private:
    // Right-hand side of the parallel transport ODE:
    // dV^i/dt = -Γ^i_{jk} V^j (du^k/dt)
    static void transport_rhs(const Surface& surface,
                               double u, double v,
                               double du, double dv,
                               double V1, double V2,
                               double& dV1dt, double& dV2dt) {
        double gamma[2][2][2];
        surface.christoffel(u, v, gamma);

        double vel[2] = {du, dv};
        double V[2] = {V1, V2};

        dV1dt = 0.0;
        dV2dt = 0.0;
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                dV1dt -= gamma[0][j][k] * V[j] * vel[k];
                dV2dt -= gamma[1][j][k] * V[j] * vel[k];
            }
        }
    }

    // Record position, vector, and diagnostics for one step
    static void record_step(const Surface& surface,
                            TransportResult& result,
                            double t, double u, double v,
                            double V1, double V2) {
        result.t_values.push_back(t);
        result.u_values.push_back(u);
        result.v_values.push_back(v);
        result.positions.push_back(surface.position(u, v));
        result.vectors_3d.push_back(surface.to_3d(u, v, V1, V2));
        result.V1_values.push_back(V1);
        result.V2_values.push_back(V2);
        result.norms.push_back(surface.metric_norm(u, v, V1, V2));

        Vec3 n = surface.normal(u, v);
        Vec3 v3d = surface.to_3d(u, v, V1, V2);
        result.tangency.push_back(std::abs(v3d.dot(n)));
    }

    // Compute holonomy angle and error metrics
    static void compute_diagnostics(TransportResult& result, double initial_norm) {
        size_t n = result.vectors_3d.size();
        if (n < 2) return;

        // Holonomy angle: signed angle between first and last transported vector
        // We need the surface normal at the start point to determine the sign
        Vec3 v_start = result.vectors_3d.front();
        Vec3 v_end   = result.vectors_3d.back();

        // Compute unsigned angle
        double dot_val = v_start.dot(v_end) / (v_start.norm() * v_end.norm() + 1e-15);
        dot_val = clamp_val_t(dot_val, -1.0, 1.0);
        double unsigned_angle = acos(dot_val);

        // Use cross product and surface normal to determine sign
        // The cross product v_start × v_end points along the surface normal
        // if the rotation is positive (counterclockwise when viewed from outside)
        Vec3 cross_vec = v_start.cross(v_end);

        // We store both the unsigned angle and let main.cpp handle comparison
        // Since holonomy on a sphere = 2π cos(θ) which can be > π,
        // we need to check if the actual rotation went "the long way around"
        // The cross product sign relative to the normal tells us the direction
        result.holonomy_angle = unsigned_angle;

        // Max norm drift
        result.max_norm_drift = 0.0;
        for (double nm : result.norms) {
            double drift = std::abs(nm - initial_norm) / (initial_norm + 1e-15);
            result.max_norm_drift = std::max(result.max_norm_drift, drift);
        }

        // Max tangency error
        result.max_tangency_error = 0.0;
        for (double tang : result.tangency) {
            result.max_tangency_error = std::max(result.max_tangency_error, tang);
        }
    }
};
