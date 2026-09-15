#pragma once
#include <array>
#include <cmath>
namespace eks {
// Independently evaluate rho for the constraint diagnostic. g is conformal.
inline double rho_scalar(const std::array<double,6>& g, double chi,
                         const std::array<double,3>& dphi, double pi,
                         double potential) {
    const double a=g[0], b=g[1], c=g[2], d=g[3], e=g[4], f=g[5];
    const double det=a*(d*f-e*e)-b*(b*f-c*e)+c*(b*e-c*d);
    const double x=dphi[0],y=dphi[1],z=dphi[2];
    const double grad=( (d*f-e*e)*x*x+(a*f-c*c)*y*y+(a*d-b*b)*z*z
        +2*(c*e-b*f)*x*y+2*(b*e-c*d)*x*z+2*(b*c-a*e)*y*z )/det;
    return 0.5*(pi*pi+chi*grad)+potential;
}
}
