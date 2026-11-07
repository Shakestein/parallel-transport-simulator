#pragma once

#include "vec3.h"
#include <string>
#include <cmath>
#include <algorithm>
#include <functional>
#include <memory>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ============================================================================
// Abstract base class for parametric surfaces r(u, v)
// ============================================================================
class Surface {
public:
    virtual ~Surface() = default;

    // Surface name for display
    virtual std::string name() const = 0;

    // Position vector r(u, v) in R^3
    virtual Vec3 position(double u, double v) const = 0;

    // Tangent vector ∂r/∂u
    virtual Vec3 r_u(double u, double v) const = 0;

    // Tangent vector ∂r/∂v
    virtual Vec3 r_v(double u, double v) const = 0;

    // Christoffel symbols of the second kind: gamma[i][j][k] = Γ^i_{jk}
    // Indices: i,j,k ∈ {0,1} where 0=u, 1=v
    virtual void christoffel(double u, double v, double gamma[2][2][2]) const = 0;

    // Parameter domain for mesh generation
    virtual void domain(double& u_min, double& u_max,
                        double& v_min, double& v_max) const = 0;

    // --- Derived methods (non-virtual) ---

    // Unit normal vector
    Vec3 normal(double u, double v) const {
        return r_u(u, v).cross(r_v(u, v)).normalized();
    }

    // First fundamental form coefficients E, F, G
    void metric(double u, double v, double& E, double& F, double& G) const {
        Vec3 ru = r_u(u, v);
        Vec3 rv = r_v(u, v);
        E = ru.dot(ru);
        F = ru.dot(rv);
        G = rv.dot(rv);
    }

    // Metric tensor as 2x2 array
    void metric_tensor(double u, double v, double g[2][2]) const {
        double E, F, G;
        metric(u, v, E, F, G);
        g[0][0] = E; g[0][1] = F;
        g[1][0] = F; g[1][1] = G;
    }

    // Convert intrinsic vector (V1, V2) to extrinsic 3D vector
    Vec3 to_3d(double u, double v, double V1, double V2) const {
        return r_u(u, v) * V1 + r_v(u, v) * V2;
    }

    // Riemannian norm: ||V||_g = sqrt(E*V1^2 + 2F*V1*V2 + G*V2^2)
    double metric_norm(double u, double v, double V1, double V2) const {
        double E, F, G;
        metric(u, v, E, F, G);
        double val = E * V1 * V1 + 2.0 * F * V1 * V2 + G * V2 * V2;
        return std::sqrt(std::max(val, 0.0));
    }
};

// ============================================================================
// Sphere of radius R
// u = θ (colatitude, 0 to π), v = φ (longitude, 0 to 2π)
// ============================================================================
class Sphere : public Surface {
public:
    double R;
    Sphere(double radius = 1.0) : R(radius) {}

    std::string name() const override { return "Sphere (R=" + std::to_string(R) + ")"; }

    void domain(double& u_min, double& u_max,
                double& v_min, double& v_max) const override {
        u_min = 0.05; u_max = M_PI - 0.05;  // avoid poles
        v_min = 0.0;  v_max = 2.0 * M_PI;
    }

    Vec3 position(double u, double v) const override {
        return {R * sin(u) * cos(v),
                R * sin(u) * sin(v),
                R * cos(u)};
    }

    Vec3 r_u(double u, double v) const override {
        return {R * cos(u) * cos(v),
                R * cos(u) * sin(v),
                -R * sin(u)};
    }

    Vec3 r_v(double u, double v) const override {
        return {-R * sin(u) * sin(v),
                R * sin(u) * cos(v),
                0.0};
    }

    void christoffel(double u, double v, double gamma[2][2][2]) const override {
        // Zero-initialize
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                for (int k = 0; k < 2; k++)
                    gamma[i][j][k] = 0.0;

        double sin_u = sin(u);
        double cos_u = cos(u);

        // Γ^θ_{φφ} = -sinθ cosθ
        gamma[0][1][1] = -sin_u * cos_u;

        // Γ^φ_{θφ} = Γ^φ_{φθ} = cotθ
        double cot_u = cos_u / (sin_u + 1e-15);
        gamma[1][0][1] = cot_u;
        gamma[1][1][0] = cot_u;
    }
};

// ============================================================================
// Torus with major radius R and minor radius r
// u = poloidal angle (around tube), v = toroidal angle (around z-axis)
// ============================================================================
class Torus : public Surface {
public:
    double R_major, r_minor;
    Torus(double R = 2.0, double r = 0.7) : R_major(R), r_minor(r) {}

    std::string name() const override {
        return "Torus (R=" + std::to_string(R_major) + ", r=" + std::to_string(r_minor) + ")";
    }

    void domain(double& u_min, double& u_max,
                double& v_min, double& v_max) const override {
        u_min = 0.0; u_max = 2.0 * M_PI;
        v_min = 0.0; v_max = 2.0 * M_PI;
    }

    Vec3 position(double u, double v) const override {
        double c = R_major + r_minor * cos(u);
        return {c * cos(v), c * sin(v), r_minor * sin(u)};
    }

    Vec3 r_u(double u, double v) const override {
        return {-r_minor * sin(u) * cos(v),
                -r_minor * sin(u) * sin(v),
                r_minor * cos(u)};
    }

    Vec3 r_v(double u, double v) const override {
        double c = R_major + r_minor * cos(u);
        return {-c * sin(v), c * cos(v), 0.0};
    }

    void christoffel(double u, double v, double gamma[2][2][2]) const override {
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                for (int k = 0; k < 2; k++)
                    gamma[i][j][k] = 0.0;

        double sin_u = sin(u);
        double cos_u = cos(u);
        double c = R_major + r_minor * cos_u;

        // Γ^u_{vv} = (R + r cosθ) sinθ / r
        gamma[0][1][1] = c * sin_u / r_minor;

        // Γ^v_{uv} = Γ^v_{vu} = -r sinθ / (R + r cosθ)
        gamma[1][0][1] = -r_minor * sin_u / (c + 1e-15);
        gamma[1][1][0] = gamma[1][0][1];
    }
};

// ============================================================================
// Hyperbolic Paraboloid (saddle surface): z = (u² - v²) / 2
// ============================================================================
class HyperbolicParaboloid : public Surface {
public:
    std::string name() const override { return "Hyperbolic Paraboloid"; }

    void domain(double& u_min, double& u_max,
                double& v_min, double& v_max) const override {
        u_min = -1.5; u_max = 1.5;
        v_min = -1.5; v_max = 1.5;
    }

    Vec3 position(double u, double v) const override {
        return {u, v, 0.5 * (u * u - v * v)};
    }

    Vec3 r_u(double u, double v) const override {
        return {1.0, 0.0, u};
    }

    Vec3 r_v(double u, double v) const override {
        return {0.0, 1.0, -v};
    }

    void christoffel(double u, double v, double gamma[2][2][2]) const override {
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                for (int k = 0; k < 2; k++)
                    gamma[i][j][k] = 0.0;

        double D = 1.0 + u * u + v * v;

        gamma[0][0][0] = u / D;     // Γ^u_{uu}
        gamma[0][1][1] = -u / D;    // Γ^u_{vv}
        gamma[1][0][0] = -v / D;    // Γ^v_{uu}
        gamma[1][1][1] = v / D;     // Γ^v_{vv}
    }
};

// ============================================================================
// Monkey Saddle: z = u³ - 3uv²
// Christoffel symbols computed via finite differences
// ============================================================================
class MonkeySaddle : public Surface {
public:
    std::string name() const override { return "Monkey Saddle"; }

    void domain(double& u_min, double& u_max,
                double& v_min, double& v_max) const override {
        u_min = -1.2; u_max = 1.2;
        v_min = -1.2; v_max = 1.2;
    }

    Vec3 position(double u, double v) const override {
        return {u, v, u * u * u - 3.0 * u * v * v};
    }

    Vec3 r_u(double u, double v) const override {
        return {1.0, 0.0, 3.0 * u * u - 3.0 * v * v};
    }

    Vec3 r_v(double u, double v) const override {
        return {0.0, 1.0, -6.0 * u * v};
    }

    void christoffel(double u, double v, double gamma[2][2][2]) const override {
        // For Monge patch z = h(u,v): Γ^i_{jk} = h_i * h_{jk} / (1 + h_u² + h_v²)
        // h = u³ - 3uv²
        // h_u = 3u² - 3v², h_v = -6uv
        // h_uu = 6u, h_uv = -6v, h_vv = -6u

        double h_u = 3.0 * u * u - 3.0 * v * v;
        double h_v = -6.0 * u * v;
        double h_uu = 6.0 * u;
        double h_uv = -6.0 * v;
        double h_vv = -6.0 * u;
        double D = 1.0 + h_u * h_u + h_v * h_v;

        // Γ^u_{jk} = h_u * h_{jk} / D
        gamma[0][0][0] = h_u * h_uu / D;
        gamma[0][0][1] = h_u * h_uv / D;
        gamma[0][1][0] = h_u * h_uv / D;
        gamma[0][1][1] = h_u * h_vv / D;

        // Γ^v_{jk} = h_v * h_{jk} / D
        gamma[1][0][0] = h_v * h_uu / D;
        gamma[1][0][1] = h_v * h_uv / D;
        gamma[1][1][0] = h_v * h_uv / D;
        gamma[1][1][1] = h_v * h_vv / D;
    }
};
