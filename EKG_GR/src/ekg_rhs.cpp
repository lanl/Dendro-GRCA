// Native BSSNCtx/ETS owns the mesh, stages, communication and AMR.
// This file supplies ONLY block derivatives, generated equations and BC/KO.
#include "ekg_native.h"
#include "ekg_derivative_specs.h"
#include "parameters.h"
#include "grDef.h"
#include "derivs.h"
#include "rhs.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>
namespace ekg {
namespace {
using Work=std::vector<std::vector<double>>;
using D1=void(*)(double* const,const double* const,double,const unsigned*,unsigned);
thread_local Work scratch;
thread_local std::array<std::vector<double>,3> ko_work;
void compute_derivatives(const double** in,unsigned offset,const double* h,
                         const unsigned* sz,unsigned bflag,Work& work) {
    const std::size_t n=static_cast<std::size_t>(sz[0])*sz[1]*sz[2];
    work.resize(generated::derivatives.size());for(auto& v:work)v.resize(n);
    const D1 d[3]={deriv_x,deriv_y,deriv_z},d2[3]={deriv_xx,deriv_yy,deriv_zz};
    for(std::size_t k=0;k<generated::derivatives.size();++k) {
        const auto p=generated::derivatives[k];
        if(p.kind==0)d[p.i](work[k].data(),in[p.field]+offset,h[p.i],sz,bflag);
        else if(p.kind==1) {
            if(p.i==p.j)d2[p.i](work[k].data(),in[p.field]+offset,h[p.i],sz,bflag);
            else d[p.j](work[k].data(),work[p.parent].data(),h[p.j],sz,bflag);
        } else {
#if EKG_ADVECTION_UPWIND
            using AD=void(*)(double* const,const double* const,double,const unsigned*,const double* const,unsigned);
            const AD ad[3]={deriv644adv_x,deriv644adv_y,deriv644adv_z};
            ad[p.i](work[k].data(),in[p.field]+offset,h[p.i],sz,in[6+p.i]+offset,bflag);
#else
            throw std::logic_error("Generated derivative manifest and compile-time advection disagree");
#endif
        }
    }
}
unsigned gradient_index(unsigned field,unsigned direction) {
    for(unsigned k=0;k<generated::derivatives.size();++k){const auto p=generated::derivatives[k];
        if(p.kind==0 && p.field==static_cast<int>(field) && p.i==static_cast<int>(direction))return k;}
    throw std::logic_error("Generated boundary gradient missing");
}
template<class F> void interior(const unsigned* sz,F&& f) {
    const unsigned pw=bssn::BSSN_PADDING_WIDTH;
    for(unsigned k=pw;k<sz[2]-pw;++k)for(unsigned j=pw;j<sz[1]-pw;++j)
        for(unsigned i=pw;i<sz[0]-pw;++i)f(i,j,k,i+sz[0]*(j+sz[1]*k));
}
bool boundary_point(unsigned i,unsigned j,unsigned k,const unsigned* s,unsigned flag) {
    const unsigned p=bssn::BSSN_PADDING_WIDTH;
    return ((flag&(1u<<OCT_DIR_LEFT)) && i==p) || ((flag&(1u<<OCT_DIR_RIGHT)) && i==s[0]-p-1)
        || ((flag&(1u<<OCT_DIR_DOWN)) && j==p) || ((flag&(1u<<OCT_DIR_UP)) && j==s[1]-p-1)
        || ((flag&(1u<<OCT_DIR_BACK)) && k==p) || ((flag&(1u<<OCT_DIR_FRONT)) && k==s[2]-p-1);
}
void boundaries_and_ko(double** out,const double** in,unsigned offset,const double* lo,
                      const double* hi,const double* h,const unsigned* sz,unsigned flag,
                      double time,const Work& work) {
    const std::size_t n=static_cast<std::size_t>(sz[0])*sz[1]*sz[2];
    for(auto& v:ko_work)v.resize(n);
    for(unsigned f=0;f<26;++f) {
        if(EKG_FREEZE_GEOMETRY && f<24)continue;
        const double sigma=f<24?bssn::KO_DISS_SIGMA:config.scalar_ko;
        if(sigma!=0) {
            ko_deriv_x(ko_work[0].data(),in[f]+offset,h[0],sz,flag);
            ko_deriv_y(ko_work[1].data(),in[f]+offset,h[1],sz,flag);
            ko_deriv_z(ko_work[2].data(),in[f]+offset,h[2],sz,flag);
            interior(sz,[&](unsigned,unsigned,unsigned,unsigned p){
                out[f][offset+p]+=sigma*(ko_work[0][p]+ko_work[1][p]+ko_work[2][p]);});
        }
    }
    if(flag && config.boundary=="sommerfeld") {
        for(unsigned f=0;f<26;++f) {
            if(EKG_FREEZE_GEOMETRY && f<24)continue;
            const double asym=(f==0||f==1||f==12||f==15||f==17)?1.0:0.0;
            const double falloff=((f>=3 && f<=5)||(f>=18 && f<=23))?2.0:1.0;
            bssn_bcs(out[f]+offset,in[f]+offset,
                work[gradient_index(f,0)].data(),work[gradient_index(f,1)].data(),
                work[gradient_index(f,2)].data(),lo,hi,falloff,asym,sz,flag);
        }
    } else if(flag && config.boundary=="analytic") {
        interior(sz,[&](unsigned i,unsigned j,unsigned k,unsigned p){
            if(!boundary_point(i,j,k,sz,flag))return;
            const auto ref=exact_solution(bssn::BSSN_ID_TYPE,{lo[0]+i*h[0],lo[1]+j*h[1],lo[2]+k*h[2]},
                                          time-bssn::BSSN_RK_TIME_BEGIN);
            for(unsigned f=0;f<26;++f)out[f][offset+p]=ref.rhs[f];
        });
    }
    if(EKG_FREEZE_GEOMETRY)
        for(unsigned f=0;f<24;++f)std::fill_n(out[f]+offset,n,0.0);
    else if(EKG_GAUGE_SYNCHRONOUS) {
        for(unsigned f:{0u,6u,7u,8u,9u,10u,11u})std::fill_n(out[f]+offset,n,0.0);
    }
}
}
void rhs_block(double** out,const double** in,unsigned offset,const double* pmin,
               const double* pmax,const unsigned* sz,unsigned bflag,double stage_time) {
    const std::size_t n=static_cast<std::size_t>(sz[0])*sz[1]*sz[2];
    double h[3];for(unsigned i=0;i<3;++i)h[i]=(pmax[i]-pmin[i])/(sz[i]-1);
    for(unsigned f=0;f<26;++f)std::fill_n(out[f]+offset,n,0.0);
    Work& work=scratch;compute_derivatives(in,offset,h,sz,bflag,work);
#include "ekg_input_aliases.inc"
#include "ekg_output_aliases.inc"
#include "ekg_derivative_aliases.inc"
    const double scalar_mu=config.mu,scalar_fa=config.fa,eta=bssn::ETA_CONST;
    constexpr double ekg_pi=3.141592653589793238462643383279502884;
    const auto lambda0=bssn::BSSN_LAMBDA[0],lambda1=bssn::BSSN_LAMBDA[1];
    const auto lambda2=bssn::BSSN_LAMBDA[2],lambda3=bssn::BSSN_LAMBDA[3];
    const auto lambda_f0=bssn::BSSN_LAMBDA_F[0],lambda_f1=bssn::BSSN_LAMBDA_F[1];
    interior(sz,[&](unsigned,unsigned,unsigned,unsigned pp){
#if EKG_FREEZE_GEOMETRY
        #include "ekg_matter_rhs.inc"
#elif EKG_FUSED_RHS
        #include "ekg_full_rhs.inc"
#else
        {
            #include "ekg_vacuum_rhs.inc"
        }
        #include "ekg_matter_rhs.inc"
        #include "ekg_apply_sources.inc"
#endif
    });
    boundaries_and_ko(out,in,offset,pmin,pmax,h,sz,bflag,stage_time,work);
}
void constraints_block(double** out,const double** in,unsigned offset,const double* pmin,
                       const double* pmax,const unsigned* sz,unsigned bflag) {
    const std::size_t n=static_cast<std::size_t>(sz[0])*sz[1]*sz[2];
    double h[3];for(unsigned i=0;i<3;++i)h[i]=(pmax[i]-pmin[i])/(sz[i]-1);
    for(unsigned f=0;f<bssn::BSSN_CONSTRAINT_NUM_VARS;++f)std::fill_n(out[f]+offset,n,0.0);
    Work& work=scratch;compute_derivatives(in,offset,h,sz,bflag,work);
#include "ekg_input_aliases.inc"
#include "ekg_derivative_aliases.inc"
    double* constraint_ham=out[0]+offset;
    double* constraint_mom0=out[1]+offset;
    double* constraint_mom1=out[2]+offset;
    double* constraint_mom2=out[3]+offset;
    const double scalar_mu=config.mu,scalar_fa=config.fa;
    constexpr double ekg_pi=3.141592653589793238462643383279502884;
    interior(sz,[&](unsigned,unsigned,unsigned,unsigned pp){
#include "ekg_constraints.inc"
    });
    // Native slots 4,5 are deliberately left zero and are NOT reported as Psi4.
    // Runtime validation disables GW extraction and excludes these from VTU.
}
}
