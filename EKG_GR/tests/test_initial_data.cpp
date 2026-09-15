#include "ekg_config.h"
#include "ekg_potential.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
namespace ekg {Config config;}
namespace {
constexpr double pi=3.141592653589793238462643383279502884;
void require(bool x,const char* m){if(!x)throw std::runtime_error(m);}
}
int main() {
    using namespace ekg;
    config.mu=.2;config.amplitude=.02;config.width=2;config.support_radius=8;
    config.id_outer_radius=24;config.id_radial_points=4096;
    prepare_initial_data(103);
    const auto& a=radial_initial_data();
    require(std::abs(a.residual)<config.id_tolerance,"Robin matching failed");
    require(a.mass>0,"Expected positive ADM mass");
    double max_res=0;
    const double h=a.h;
    for(unsigned i=4;i+2<a.psi.size();++i) {
        const double r=i*h,p=a.psi[i];
        const double dp=(a.psi[i-2]-8*a.psi[i-1]+8*a.psi[i+1]-a.psi[i+2])/(12*h);
        const double dd=(-a.psi[i+2]+16*a.psi[i+1]-30*p+16*a.psi[i-1]-a.psi[i-2])/(12*h*h);
        const auto ph=radial_profile(r);
        const double H=-8*std::pow(p,-5)*(dd+2*dp/r)-8*pi*std::pow(p,-4)*ph[1]*ph[1]
                       -16*pi*potential(ph[0],config.mu,config.fa);
        max_res=std::max(max_res,std::abs(H));
    }
    require(max_res<2e-6,"Radial finite-difference Hamiltonian residual too large");
    require(radial_profile(config.support_radius)[0]==0,"Profile not compact");
    config.mu=0;config.pi0=.05;
    for(double t:{0.,.2,1.}) {
        const auto e=exact_solution(102,{0,0,0},t);
        require(std::abs(2*e.state[2]*e.state[2]/3-8*pi*e.state[25]*e.state[25])<1e-13,"FLRW Hamiltonian mismatch");
        require(std::abs(e.rhs[25]-e.state[2]*e.state[25])<1e-13,"FLRW Pi evolution mismatch");
    }
    std::cout<<"PASS: nonlinear radial Hamiltonian shooting and compact scalar data\n"
             <<"ADM mass="<<a.mass<<", Robin residual="<<a.residual
             <<", radial FD max|H|="<<max_res<<"\n"
             <<"PASS: exact coupled FLRW constraints/evolution\n";
}
