#pragma once
#include <array>
#include <cmath>
#include <stdexcept>
namespace eks {
// All coordinates in this header are PHYSICAL, not octree coordinates.
struct ExactState { double phi, pi, pi_t; };
inline ExactState gaussian_exact(double r, double t, double A, double s) {
    if (!(s > 0.0)) throw std::invalid_argument("Gaussian width must be positive");
    const double s2=s*s;
    if (std::abs(r) < 1.0e-6*s) {
        // Analytic r -> 0 limit. Avoid dividing two cancelling differences.
        const double u=t*t/s2, e=std::exp(-0.5*u);
        return {A*(1.0-u)*e,
                A*t/s2*(3.0-u)*e,
                A/s2*(3.0-6.0*u+u*u)*e};
    }
    const double a=r-t, b=r+t;
    const double ea=std::exp(-a*a/(2.0*s2));
    const double eb=std::exp(-b*b/(2.0*s2));
    return {A*(a*ea+b*eb)/(2.0*r),
            A*((1.0-a*a/s2)*ea-(1.0-b*b/s2)*eb)/(2.0*r),
            -A*((a*a*a/(s2*s2)-3.0*a/s2)*ea+
                (b*b*b/(s2*s2)-3.0*b/s2)*eb)/(2.0*r)};
}
inline ExactState plane_exact(const std::array<double,3>& x, double t,
                              double A, double mu,
                              const std::array<double,3>& k) {
    const double w=std::sqrt(mu*mu+k[0]*k[0]+k[1]*k[1]+k[2]*k[2]);
    const double q=k[0]*x[0]+k[1]*x[1]+k[2]*x[2]-w*t;
    return {A*std::cos(q), -A*w*std::sin(q), A*w*w*std::cos(q)};
}
}
