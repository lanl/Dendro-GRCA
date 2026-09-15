#include "ekg_config.h"
#include "ekg_potential.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace ekg {
namespace {
constexpr double pi=3.141592653589793238462643383279502884;
RadialData radial;
State flat() { State u{}; u[0]=u[1]=u[12]=u[15]=u[17]=1.0;return u; }
double radius(const std::array<double,3>& x) {
    double r2=0;for(unsigned i=0;i<3;++i)r2+=(x[i]-config.center[i])*(x[i]-config.center[i]);
    return std::sqrt(r2);
}
using Y=std::array<double,2>;
Y ode(double r,const Y& y) {
    const auto f=radial_profile(r);
    const double source=-pi*y[0]*f[1]*f[1]-2*pi*std::pow(y[0],5)*potential(f[0],config.mu,config.fa);
    return {y[1],r==0?source/3:source-2*y[1]/r};
}
Y step(double r,const Y& y,double h) {
    auto plus=[](const Y& a,const Y& b,double f){return Y{a[0]+f*b[0],a[1]+f*b[1]};};
    const auto k1=ode(r,y),k2=ode(r+h/2,plus(y,k1,h/2));
    const auto k3=ode(r+h/2,plus(y,k2,h/2)),k4=ode(r+h,plus(y,k3,h));
    return {y[0]+h*(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6,
            y[1]+h*(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6};
}
RadialData integrate(double central) {
    RadialData a;a.outer=config.id_outer_radius;a.h=a.outer/config.id_radial_points;
    const unsigned n=config.id_radial_points;
    a.psi.resize(n+1);a.dpsi.resize(n+1);
    a.psi[0]=central;a.dpsi[0]=0;
    const double src0=ode(0,{central,0})[1]*3;
    Y y{central+src0*a.h*a.h/6,src0*a.h/3};
    for(unsigned k=1;k<=n;++k) {
        if(!std::isfinite(y[0])||!std::isfinite(y[1])||y[0]<=0)
            throw std::runtime_error("Radial Hamiltonian shoot lost positive conformal factor");
        a.psi[k]=y[0];a.dpsi[k]=y[1];
        if(k<n)y=step(k*a.h,y,a.h);
    }
    a.residual=a.psi.back()+a.outer*a.dpsi.back()-1;
    a.mass=-2*a.outer*a.outer*a.dpsi.back();
    return a;
}
}
std::array<double,2> radial_profile(double r) {
    const double R=config.support_radius;
    if(r>=R)return {0,0};
    const double den=R*R-r*r;
    const double p=config.amplitude*std::exp(-r*r/(2*config.width*config.width)-r*r/den);
    if(p==0)return {0,0};
    return {p,p*(-r/(config.width*config.width)-2*r*R*R/(den*den))};
}
double RadialData::value(double r) const {
    if(psi.empty())throw std::logic_error("Radial initial data not prepared");
    if(r>=outer)return 1+mass/(2*r);
    const auto i=std::min<std::size_t>(static_cast<std::size_t>(r/h),psi.size()-2);
    const double t=(r-i*h)/h,t2=t*t,t3=t2*t;
    return (2*t3-3*t2+1)*psi[i]+(t3-2*t2+t)*h*dpsi[i]
           +(-2*t3+3*t2)*psi[i+1]+(t3-t2)*h*dpsi[i+1];
}
const RadialData& radial_initial_data(){return radial;}
void prepare_initial_data(unsigned id) {
    if(id!=103)return;
    if(config.id_outer_radius<=config.support_radius || config.id_radial_points<128)
        throw std::runtime_error("ID103 needs id_outer_radius>support_radius and >=128 radial points");
    auto lo_solution=integrate(1.0);
    if(std::abs(lo_solution.residual)<config.id_tolerance){radial=std::move(lo_solution);return;}
    double lo=1.0,hi=lo,fl=lo_solution.residual;
    bool bracket=false;
    for(unsigned k=1;k<=200;++k) {
        hi=1.0+0.025*k;
        const auto shot=integrate(hi);
        if(fl*shot.residual<=0){bracket=true;break;}
    }
    if(!bracket)throw std::runtime_error("No weak positive branch bracketed for ID103; reduce scalar amplitude or solve elliptic data with a more general solver");
    for(unsigned n=0;n<80;++n) {
        const double mid=(lo+hi)/2;auto shot=integrate(mid);
        if(std::abs(shot.residual)<config.id_tolerance){radial=std::move(shot);return;}
        if(fl*shot.residual<=0)hi=mid;else {lo=mid;fl=shot.residual;}
    }
    throw std::runtime_error("Radial Hamiltonian shooting tolerance not reached");
}
Exact exact_solution(unsigned id,const std::array<double,3>& x,double t) {
    Exact z;z.state=flat();
    if(id==100) {
        const double r=radius(x),s2=config.width*config.width,A=config.amplitude;
        auto F=[&](double q){return q*std::exp(-q*q/(2*s2));};
        auto D=[&](double q){return (1-q*q/s2)*std::exp(-q*q/(2*s2));};
        auto D2=[&](double q){return (q*q*q/(s2*s2)-3*q/s2)*std::exp(-q*q/(2*s2));};
        if(r<1e-7*config.width) {
            const double e=std::exp(-t*t/(2*s2));
            z.state[24]=A*(1-t*t/s2)*e;
            z.state[25]=A*(3*t/s2-t*t*t/(s2*s2))*e;
            z.rhs[25]=A*(3/s2-6*t*t/(s2*s2)+std::pow(t,4)/(s2*s2*s2))*e;
        } else {
            z.state[24]=A*(F(r-t)+F(r+t))/(2*r);
            z.state[25]=A*(D(r-t)-D(r+t))/(2*r);
            z.rhs[25]=-A*(D2(r-t)+D2(r+t))/(2*r);
        }
        z.rhs[24]=-z.state[25];
    } else if(id==101) {
        double k2=0,theta=0;
        for(unsigned i=0;i<3;++i){k2+=config.wavevector[i]*config.wavevector[i];theta+=config.wavevector[i]*(x[i]-config.center[i]);}
        const double w=std::sqrt(k2+config.mu*config.mu);theta-=w*t;
        z.state[24]=config.amplitude*std::cos(theta);
        z.state[25]=-config.amplitude*w*std::sin(theta);
        z.rhs[24]=-z.state[25];z.rhs[25]=w*w*z.state[24];
    } else if(id==102) {
        if(config.pi0==0)throw std::runtime_error("FLRW requires nonzero pi0");
        const double t0=1/(std::sqrt(12*pi)*std::abs(config.pi0)),u=1+t/t0;
        if(u<=0)throw std::runtime_error("FLRW time outside expanding branch");
        z.state[1]=std::pow(u,-2.0/3);z.state[2]=-1/(t0+t);
        z.state[24]=config.phi0-config.pi0*t0*std::log(u);z.state[25]=config.pi0/u;
        z.rhs[1]=-2*z.state[1]/(3*(t0+t));z.rhs[2]=1/((t0+t)*(t0+t));
        z.rhs[24]=-z.state[25];z.rhs[25]=-config.pi0/(t0*u*u);
    } else throw std::runtime_error("No analytic time-dependent reference for this ID");
    return z;
}
State initial_state(unsigned id,const std::array<double,3>& x) {
    if(id>=100 && id<=102)return exact_solution(id,x,0).state;
    if(id==103) {
        State u=flat();const double r=radius(x),psi=radial.value(r);
        u[1]=1/std::pow(psi,4);u[24]=radial_profile(r)[0];u[25]=0;
        return u;
    }
    throw std::runtime_error("Not an EKG-specific initial data ID");
}
} // namespace ekg
