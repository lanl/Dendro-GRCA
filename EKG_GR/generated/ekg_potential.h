#pragma once
#include <cmath>
namespace ekg {
inline double potential(double phi,double mu,double fa) {
    return (1.0/2.0)*pow(mu, 2)*pow(phi, 2);
}
inline double potential_prime(double phi,double mu,double fa) {
    return pow(mu, 2)*phi;
}
}
