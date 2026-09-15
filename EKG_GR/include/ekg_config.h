#pragma once
#include <array>
#include <string>
#include <vector>
namespace ekg {
inline constexpr unsigned field_count=26;
inline constexpr unsigned phi_index=24, pi_index=25;
struct Config {
    double mu=0.2, fa=1.0;
    double amplitude=0.02, width=2.0, support_radius=8.0;
    double phi0=0.0, pi0=0.05;
    std::array<double,3> center{0,0,0}, wavevector{0.5,0.25,0};
    std::string boundary="sommerfeld";
    double scalar_ko=0.0;
    // Fixed refinement normalization, never the instantaneous amplitude.
    double refine_scale_phi=1.0, refine_scale_pi=1.0;
    double id_outer_radius=24.0, id_tolerance=1e-12;
    unsigned id_radial_points=4096;
    bool allow_constraint_violating_id=false;
};
extern Config config;
using State=std::array<double,field_count>;
struct Exact { State state{}, rhs{}; };
Exact exact_solution(unsigned id,const std::array<double,3>& x,double elapsed);
State initial_state(unsigned id,const std::array<double,3>& x);
void prepare_initial_data(unsigned id);
// Compact C-infinity radial field and its radial derivative.
std::array<double,2> radial_profile(double r);
struct RadialData {
    double h=0, outer=0, mass=0, residual=0;
    std::vector<double> psi,dpsi;
    double value(double r) const;
};
const RadialData& radial_initial_data();
} // namespace ekg
