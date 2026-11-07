#pragma once

#include <cmath>
#include <string>
#include <functional>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Helper: clamp value to [lo, hi] — compatible with older compilers
inline double clamp_val(double x, double lo, double hi) {
    return (x < lo) ? lo : (x > hi) ? hi : x;
}

// ============================================================================
// Abstract base class for parametric curves in (u, v) parameter space
// ============================================================================
class Curve {
public:
    virtual ~Curve() = default;

    virtual std::string name() const = 0;

    // Evaluate curve position and velocity at parameter t
    // Returns: u(t), v(t), du/dt, dv/dt
    virtual void evaluate(double t, double& u, double& v,
                          double& du, double& dv) const = 0;

    virtual double t_start() const = 0;
    virtual double t_end() const = 0;
};

// ============================================================================
// Latitude Circle on a sphere: θ = θ₀ (constant), φ goes from 0 to 2π
// Classic holonomy demonstration
// ============================================================================
class LatitudeCircle : public Curve {
public:
    double theta0;  // colatitude

    LatitudeCircle(double theta = M_PI / 4.0) : theta0(theta) {}

    std::string name() const override {
        return "Latitude Circle (theta=" + std::to_string(theta0) + ")";
    }

    void evaluate(double t, double& u, double& v,
                  double& du, double& dv) const override {
        u = theta0;
        v = t;
        du = 0.0;
        dv = 1.0;
    }

    double t_start() const override { return 0.0; }
    double t_end() const override { return 2.0 * M_PI; }
};

// ============================================================================
// Great Circle (geodesic) on a sphere
// Parametrized as a tilted great circle
// ============================================================================
class GreatCircle : public Curve {
public:
    double tilt;  // tilt angle from equator

    GreatCircle(double tilt_angle = M_PI / 4.0) : tilt(tilt_angle) {}

    std::string name() const override {
        return "Great Circle (tilt=" + std::to_string(tilt) + ")";
    }

    void evaluate(double t, double& u, double& v,
                  double& du, double& dv) const override {
        // Great circle: point on sphere = (sin(t), cos(t)*sin(tilt), cos(t)*cos(tilt))
        // Convert to spherical: θ = acos(z), φ = atan2(y, x)
        double x = sin(t);
        double y = cos(t) * sin(tilt);
        double z = cos(t) * cos(tilt);

        u = acos(clamp_val(z, -1.0, 1.0));  // θ
        v = atan2(y, x);                       // φ
        if (v < 0) v += 2.0 * M_PI;

        // Numerical derivatives (central difference)
        double dt = 1e-7;
        double x1 = sin(t + dt), y1 = cos(t + dt) * sin(tilt), z1 = cos(t + dt) * cos(tilt);
        double x0 = sin(t - dt), y0 = cos(t - dt) * sin(tilt), z0 = cos(t - dt) * cos(tilt);

        double u1 = acos(clamp_val(z1, -1.0, 1.0));
        double v1 = atan2(y1, x1); if (v1 < 0) v1 += 2.0 * M_PI;
        double u0 = acos(clamp_val(z0, -1.0, 1.0));
        double v0 = atan2(y0, x0); if (v0 < 0) v0 += 2.0 * M_PI;

        du = (u1 - u0) / (2.0 * dt);

        // Handle φ wraparound
        double dv_raw = v1 - v0;
        if (dv_raw > M_PI) dv_raw -= 2.0 * M_PI;
        if (dv_raw < -M_PI) dv_raw += 2.0 * M_PI;
        dv = dv_raw / (2.0 * dt);
    }

    double t_start() const override { return 0.0; }
    double t_end() const override { return 2.0 * M_PI; }
};

// ============================================================================
// Spiral curve on a surface
// u(t) = u0 + a*t, v(t) = omega*t
// ============================================================================
class SpiralCurve : public Curve {
public:
    double u0, a, omega, T;

    SpiralCurve(double u_start = 0.3, double slope = 0.0,
                double angular_speed = 1.0, double duration = 2.0 * M_PI)
        : u0(u_start), a(slope), omega(angular_speed), T(duration) {}

    std::string name() const override { return "Spiral Curve"; }

    void evaluate(double t, double& u, double& v,
                  double& du, double& dv) const override {
        u = u0 + a * t;
        v = omega * t;
        du = a;
        dv = omega;
    }

    double t_start() const override { return 0.0; }
    double t_end() const override { return T; }
};

// ============================================================================
// Circular path in (u, v) parameter space
// u(t) = u_c + r*cos(t), v(t) = v_c + r*sin(t)
// ============================================================================
class CircularPath : public Curve {
public:
    double u_c, v_c, radius;

    CircularPath(double center_u = 0.0, double center_v = 0.0, double r = 0.8)
        : u_c(center_u), v_c(center_v), radius(r) {}

    std::string name() const override { return "Circular Path in Parameter Space"; }

    void evaluate(double t, double& u, double& v,
                  double& du, double& dv) const override {
        u = u_c + radius * cos(t);
        v = v_c + radius * sin(t);
        du = -radius * sin(t);
        dv = radius * cos(t);
    }

    double t_start() const override { return 0.0; }
    double t_end() const override { return 2.0 * M_PI; }
};

// ============================================================================
// Torus Knot: u(t) = p*t, v(t) = q*t
// Wraps p times poloidally and q times toroidally
// ============================================================================
class TorusKnot : public Curve {
public:
    int p, q;

    TorusKnot(int p_wraps = 1, int q_wraps = 3) : p(p_wraps), q(q_wraps) {}

    std::string name() const override {
        return "Torus Knot (" + std::to_string(p) + "," + std::to_string(q) + ")";
    }

    void evaluate(double t, double& u, double& v,
                  double& du, double& dv) const override {
        u = p * t;
        v = q * t;
        du = static_cast<double>(p);
        dv = static_cast<double>(q);
    }

    double t_start() const override { return 0.0; }
    double t_end() const override { return 2.0 * M_PI; }
};

// ============================================================================
// Toroidal Loop: constant u (poloidal angle), v goes 0 to 2π
// Analogous to a latitude circle but on a torus
// ============================================================================
class ToroidalLoop : public Curve {
public:
    double u0;

    ToroidalLoop(double poloidal_angle = 0.0) : u0(poloidal_angle) {}

    std::string name() const override {
        return "Toroidal Loop (u=" + std::to_string(u0) + ")";
    }

    void evaluate(double t, double& u, double& v,
                  double& du, double& dv) const override {
        u = u0;
        v = t;
        du = 0.0;
        dv = 1.0;
    }

    double t_start() const override { return 0.0; }
    double t_end() const override { return 2.0 * M_PI; }
};
