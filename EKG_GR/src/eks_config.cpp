#include "eks_config.h"
#include "eks_exact.h"
#include "parameters.h"
#include "grDef.h"
#include "eks_generated_metadata.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <vector>
namespace eks {
Config config;
void read_config(const toml::value& root) {
    if (!root.contains("EKS")) throw std::runtime_error("Missing [EKS] table");
    const auto& t=root.at("EKS");
    const std::set<std::string> allowed={"mode","boundary","mu","fa","amplitude",
      "width","ko_sigma","center","wavevector","allow_constraint_violating_id"};
    for (const auto& item:t.as_table())
        if (!allowed.count(item.first)) throw std::runtime_error("Unknown EKS key: "+item.first);
    config.mode=toml::find_or<std::string>(t,"mode","frozen");
    config.boundary=toml::find_or<std::string>(t,"boundary","analytic");
    config.mu=toml::find_or<double>(t,"mu",0.0);
    config.fa=toml::find_or<double>(t,"fa",1.0);
    config.amplitude=toml::find_or<double>(t,"amplitude",1.0e-3);
    config.width=toml::find_or<double>(t,"width",2.0);
    config.ko_sigma=toml::find_or<double>(t,"ko_sigma",0.0);
    config.allow_constraint_violating_id=toml::find_or<bool>(t,"allow_constraint_violating_id",false);
    for (auto name:{"center","wavevector"}) {
        if (!t.contains(name)) continue;
        auto a=toml::find<std::vector<double>>(t,name);
        if (a.size()!=3) throw std::runtime_error(std::string(name)+" requires three floats");
        auto& target=std::string(name)=="center"?config.center:config.wavevector;
        std::copy(a.begin(),a.end(),target.begin());
    }
    for (double a:{config.mu,config.fa,config.amplitude,config.width,config.ko_sigma})
        if (!std::isfinite(a)) throw std::runtime_error("Non-finite EKS parameter");
    for (double a:config.center) if (!std::isfinite(a)) throw std::runtime_error("Non-finite center");
    for (double a:config.wavevector) if (!std::isfinite(a)) throw std::runtime_error("Non-finite wavevector");
    if(config.mu<0 || config.fa<=0 || config.width<=0 || config.ko_sigma<0)
        throw std::runtime_error("Require mu>=0, fa>0, width>0, ko_sigma>=0");
    if(config.mode!="frozen" && config.mode!="vacuum_background" && config.mode!="coupled")
        throw std::runtime_error("EKS mode must be frozen, vacuum_background, or coupled");
    if(config.boundary!="analytic" && config.boundary!="sommerfeld")
        throw std::runtime_error("EKS boundary must be analytic or sommerfeld");
}
toml::value config_as_toml() {
    toml::value t=toml::table{};
    t["mode"]=config.mode; t["boundary"]=config.boundary;
    t["mu"]=config.mu; t["fa"]=config.fa; t["amplitude"]=config.amplitude;
    t["width"]=config.width; t["ko_sigma"]=config.ko_sigma;
    t["center"]=std::vector<double>(config.center.begin(),config.center.end());
    t["wavevector"]=std::vector<double>(config.wavevector.begin(),config.wavevector.end());
    t["allow_constraint_violating_id"]=config.allow_constraint_violating_id;
    return t;
}
void validate_runtime(unsigned int ts_mode) {
    static_assert(bssn::BSSN_NUM_VARS==26,"EKS requires 26 evolution fields");
    static_assert(bssn::U_SCALARPHI==24 && bssn::U_SCALARPI==25,"EKS layout mismatch");
    if(ts_mode!=1 || bssn::BSSN_RK_TYPE!=1)
        throw std::runtime_error("Starter supports UTS ts_mode=1 and RK4 only");
    if(bssn::BSSN_ELE_ORDER!=6)
        throw std::runtime_error("Starter CMake selects sixth-order derivatives: set BSSN_ELE_ORDER=6");
    if(!bssn::BSSN_ASYNC_COMM_K || 26%bssn::BSSN_ASYNC_COMM_K)
        throw std::runtime_error("BSSN_ASYNC_COMM_K must divide 26; use 1 or 2");
    if(bssn::BSSN_RESTORE_SOLVER)
        throw std::runtime_error("Starter restart is disabled until EKS metadata/schema checks are added");
    if(bssn::BSSN_ID_TYPE!=100 && bssn::BSSN_ID_TYPE!=101)
        throw std::runtime_error("Starter ships only ID 100 Gaussian and 101 plane wave; add a vetted BH ID next");
    if(config.boundary=="analytic" && config.mode!="frozen")
        throw std::runtime_error("Exact flat-space boundary requires EKS mode=frozen");
    if(bssn::BSSN_ID_TYPE==100 && config.mu!=0.0)
        throw std::runtime_error("ID100 is the exact MASSLESS Gaussian test: mu must be zero");
    if(config.boundary=="analytic" && std::string(EKS_GENERATED_POTENTIAL)!="massive")
        throw std::runtime_error("Exact tests require the free massive/massless generator");
    if(config.mode=="coupled" && !config.allow_constraint_violating_id)
        throw std::runtime_error("Flat+scalar ID violates Einstein constraints. Supply solved ID or explicitly acknowledge a DEVELOPMENT-only coupled run");
    if(config.mode=="coupled" && (AEH::AEH_SOLVER_FREQ || bssn::BSSN_GW_EXTRACT_FREQ))
        throw std::runtime_error("Starter coupled development runs do not validate horizon/GW extraction");
    if(config.mode=="coupled")
        std::cerr << "WARNING: constraint-violating development ID; NOT production Einstein-scalar data\n";
}
void initial_data_physical(double x,double y,double z,double* u) {
    std::fill_n(u,bssn::BSSN_NUM_VARS,0.0);
    u[bssn::U_ALPHA]=u[bssn::U_CHI]=1.0;
    u[bssn::U_SYMGT0]=u[bssn::U_SYMGT3]=u[bssn::U_SYMGT5]=1.0;
    const std::array<double,3> p={x-config.center[0],y-config.center[1],z-config.center[2]};
    ExactState q;
    if(bssn::BSSN_ID_TYPE==100) {
        const double r=std::sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
        q=gaussian_exact(r,bssn::BSSN_RK_TIME_BEGIN,config.amplitude,config.width);
    } else {
        q=plane_exact(p,bssn::BSSN_RK_TIME_BEGIN,config.amplitude,config.mu,config.wavevector);
    }
    u[bssn::U_SCALARPHI]=q.phi;
    u[bssn::U_SCALARPI]=q.pi;
}
}
