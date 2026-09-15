#include "ekg_build_config.h"
#include "ekg_derivative_specs.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>
using State=std::array<double,26>;
using Work=std::vector<std::vector<double>>;
State evaluate(State u,const Work& work,bool fused) {
    const unsigned offset=0,pp=0;
    std::array<const double*,26> in_storage;
    std::array<double*,26> out_storage;
    State result{};
    for(unsigned i=0;i<26;++i){in_storage[i]=&u[i];out_storage[i]=&result[i];}
    const double** in=in_storage.data();double** out=out_storage.data();
#include "ekg_input_aliases.inc"
#include "ekg_output_aliases.inc"
#include "ekg_derivative_aliases.inc"
    const double scalar_mu=0.2,scalar_fa=0.7,eta=2;
    const double lambda0=1,lambda1=1,lambda2=1,lambda3=1,lambda_f0=1,lambda_f1=0;
    constexpr double ekg_pi=3.141592653589793238462643383279502884;
    if(fused) {
#include "ekg_full_rhs.inc"
    } else {
        {
#include "ekg_vacuum_rhs.inc"
        }
#include "ekg_matter_rhs.inc"
#include "ekg_apply_sources.inc"
    }
    return result;
}
std::array<double,4> evaluate_constraints(State u,const Work& work) {
    const unsigned offset=0,pp=0;
    std::array<const double*,26> storage;
    for(unsigned i=0;i<26;++i)storage[i]=&u[i];
    const double** in=storage.data();
#include "ekg_input_aliases.inc"
#include "ekg_derivative_aliases.inc"
    std::array<double,4> result{};
    double* constraint_ham=&result[0];double* constraint_mom0=&result[1];
    double* constraint_mom1=&result[2];double* constraint_mom2=&result[3];
    const double scalar_mu=.2,scalar_fa=.7;
    constexpr double ekg_pi=3.141592653589793238462643383279502884;
#include "ekg_constraints.inc"
    return result;
}
void require(bool pass,const char* msg){if(!pass)throw std::runtime_error(msg);}
int main() {
    std::mt19937_64 rng(89234);std::uniform_real_distribution<double> ran(-.05,.05);
    Work w(ekg::generated::derivatives.size(),std::vector<double>(1));
    double worst=0;
    for(unsigned sample=0;sample<30;++sample) {
        State u{};for(double& v:u)v=ran(rng);
        u[0]=.9;u[1]=.8;u[12]+=1;u[15]+=1;u[17]+=1;
        for(auto& v:w)v[0]=ran(rng);
        const auto full=evaluate(u,w,true),split=evaluate(u,w,false);
        for(unsigned i=0;i<26;++i) {
            require(std::isfinite(full[i])&&std::isfinite(split[i]),"nonfinite generated kernel");
            worst=std::max(worst,std::abs(full[i]-split[i]));
            require(std::abs(full[i]-split[i])<1e-10*(1+std::abs(full[i])),"fused and split kernels disagree");
        }
    }
    State u{};u[0]=u[1]=u[12]=u[15]=u[17]=1;u[24]=.08;u[25]=-.03;
    for(auto& v:w)v[0]=0;
    const double k[3]={.2,.1,.3};
    for(unsigned n=0;n<w.size();++n){auto p=ekg::generated::derivatives[n];
        if(p.kind==1 && p.field==24)w[n][0]=-k[p.i]*k[p.j]*u[24];}
    const auto r=evaluate(u,w,true);
    require(std::abs(r[24]+u[25])<1e-13,"Pi=-dt(phi) convention broken");
    const double Vp=EKG_POTENTIAL_AXION?(.2*.2*.7*std::sin(u[24]/.7)):(.2*.2*u[24]);
    require(std::abs(r[25]-(.2*.2+.1*.1+.3*.3)*u[24]-Vp)<1e-12,"flat scalar Laplacian/potential sign broken");
    // Independent conformally-flat constraint calculation, nonzero scalar momentum.
    for(auto& v:w)v[0]=0;
    u[1]=.9;const double dc[3]={.02,.01,-.01},dp[3]={.03,-.02,.01};
    for(unsigned n=0;n<w.size();++n){auto p=ekg::generated::derivatives[n];
        if(p.kind==0 && p.field==1)w[n][0]=dc[p.i];
        if(p.kind==0 && p.field==24)w[n][0]=dp[p.i];
        if(p.kind==1 && p.field==1 && p.i==p.j)w[n][0]=.001;}
    const auto c=evaluate_constraints(u,w);
    constexpr double pi=3.141592653589793238462643383279502884;
    const double V=EKG_POTENTIAL_AXION ? 2*.2*.2*.7*.7*std::pow(std::sin(u[24]/(2*.7)),2)
                                          : .5*.2*.2*u[24]*u[24];
    double dc2=0,dp2=0;for(unsigned i=0;i<3;++i){dc2+=dc[i]*dc[i];dp2+=dp[i]*dp[i];}
    const double expectedH=2*.003-2.5*dc2/u[1]
        -(EKG_BACKREACTION ? (8*pi*(u[25]*u[25]+u[1]*dp2)+16*pi*V):0);
    require(std::abs(c[0]-expectedH)<1e-12,"Generated Hamiltonian constraint mismatch");
    for(unsigned i=0;i<3;++i)
        require(std::abs(c[1+i]+(EKG_BACKREACTION?8*pi*u[25]*dp[i]:0))<1e-12,"Generated momentum constraint mismatch");
    std::cout<<"PASS: independently checked conformally-flat H and all M_i\n";
    std::cout<<"PASS: 30 general-jet fused/split comparisons, 26 outputs each\n"
             <<"Maximum absolute discrepancy: "<<worst<<"\n"
             <<"PASS: generated flat-space scalar signs and potential\n";
}
