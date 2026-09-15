#include "ekg_native.h"
#include "parameters.h"
#include "grDef.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <set>
#include <stdexcept>
namespace ekg {
Config config;
void read_config(const toml::value& root) {
    if(!root.contains("EKG"))throw std::runtime_error("Missing [EKG] table");
    const auto& t=root.at("EKG");
    const std::set<std::string> allowed={"mu","fa","amplitude","width","support_radius",
        "phi0","pi0","center","wavevector","boundary","scalar_ko","refine_scale_phi",
        "refine_scale_pi","id_outer_radius","id_radial_points","id_tolerance",
        "allow_constraint_violating_id"};
    for(const auto& kv:t.as_table())if(!allowed.count(kv.first))
        throw std::runtime_error("Unknown EKG option: "+kv.first+"; coupling/gauge/potential are CMake options");
#define EKG_READ_DOUBLE(name) config.name=toml::find_or<double>(t,#name,config.name)
    EKG_READ_DOUBLE(mu); EKG_READ_DOUBLE(fa); EKG_READ_DOUBLE(amplitude);
    EKG_READ_DOUBLE(width); EKG_READ_DOUBLE(support_radius); EKG_READ_DOUBLE(phi0);
    EKG_READ_DOUBLE(pi0); EKG_READ_DOUBLE(scalar_ko); EKG_READ_DOUBLE(refine_scale_phi);
    EKG_READ_DOUBLE(refine_scale_pi); EKG_READ_DOUBLE(id_outer_radius); EKG_READ_DOUBLE(id_tolerance);
#undef EKG_READ_DOUBLE
    config.boundary=toml::find_or<std::string>(t,"boundary",config.boundary);
    config.id_radial_points=toml::find_or<unsigned>(t,"id_radial_points",config.id_radial_points);
    config.allow_constraint_violating_id=toml::find_or<bool>(t,"allow_constraint_violating_id",false);
    for(const auto* key:{"center","wavevector"})if(t.contains(key)) {
        auto v=toml::find<std::vector<double>>(t,key);
        if(v.size()!=3)throw std::runtime_error(std::string(key)+" must have three elements");
        auto& target=std::string(key)=="center"?config.center:config.wavevector;
        std::copy(v.begin(),v.end(),target.begin());
    }
    for(double v:{config.mu,config.fa,config.amplitude,config.width,config.support_radius,
        config.phi0,config.pi0,config.scalar_ko,config.refine_scale_phi,config.refine_scale_pi,
        config.id_outer_radius,config.id_tolerance})
        if(!std::isfinite(v))throw std::runtime_error("Non-finite EKG parameter");
    for(const auto& v:config.center)if(!std::isfinite(v))throw std::runtime_error("Non-finite EKG center");
    for(const auto& v:config.wavevector)if(!std::isfinite(v))throw std::runtime_error("Non-finite wavevector");
    if(config.mu<0 || config.fa<=0 || config.width<=0 || config.support_radius<=0 ||
       config.scalar_ko<0 || config.refine_scale_phi<=0 || config.refine_scale_pi<=0 || config.id_tolerance<=0)
        throw std::runtime_error("Invalid scalar/refinement/initial-data scale");
    if(config.boundary!="analytic" && config.boundary!="sommerfeld")
        throw std::runtime_error("EKG boundary must be analytic or sommerfeld");
    for(const auto* key:{"BSSN_VTU_FILE_PREFIX","BSSN_CHKPT_FILE_PREFIX",
                         "BSSN_PROFILE_FILE_PREFIX","DENDRO_LOG_FILE"}) {
        if(root.contains(key)) {
            const auto p=std::filesystem::path(toml::find<std::string>(root,key)).parent_path();
            if(!p.empty())std::filesystem::create_directories(p);
        }
    }
}
toml::value config_as_toml() {
    toml::value t=toml::table{};
#define EKG_DUMP(name) t[#name]=config.name
    EKG_DUMP(mu);EKG_DUMP(fa);EKG_DUMP(amplitude);EKG_DUMP(width);EKG_DUMP(support_radius);
    EKG_DUMP(phi0);EKG_DUMP(pi0);EKG_DUMP(scalar_ko);EKG_DUMP(refine_scale_phi);EKG_DUMP(refine_scale_pi);
    EKG_DUMP(id_outer_radius);EKG_DUMP(id_radial_points);EKG_DUMP(id_tolerance);
    EKG_DUMP(boundary);EKG_DUMP(allow_constraint_violating_id);
#undef EKG_DUMP
    t["center"]=std::vector<double>(config.center.begin(),config.center.end());
    t["wavevector"]=std::vector<double>(config.wavevector.begin(),config.wavevector.end());
    return t;
}
void validate_runtime(unsigned ts_mode,MPI_Comm comm) {
#include "ekg_field_asserts.inc"
    static_assert(EKG_BACKREACTION==EKG_GENERATED_BACKREACTION,"Stale coupling kernel");
    static_assert(EKG_GAUGE_SYNCHRONOUS==EKG_GENERATED_GAUGE_SYNCHRONOUS,"Stale gauge kernel");
    static_assert(EKG_POTENTIAL_AXION==EKG_GENERATED_POTENTIAL_AXION,"Stale potential kernel");
    static_assert(EKG_ADVECTION_UPWIND==EKG_GENERATED_ADVECTION_UPWIND,"Stale derivative kernel");
#if EKG_GENERATED_REFERENCE_EMITTER
#error "ekgSolver must use DendroSym-generated kernels, not the reference test emitter"
#endif
    if(ts_mode!=1)throw std::runtime_error("This integration targets the native ETS (ts_mode=1); LTS needs a separate validation campaign");
    if(bssn::BSSN_ELE_ORDER!=6)throw std::runtime_error("This module currently selects native sixth-order derivatives; use ELE_ORDER=6");
    if(!bssn::BSSN_ASYNC_COMM_K || 26%bssn::BSSN_ASYNC_COMM_K)
        throw std::runtime_error("BSSN_ASYNC_COMM_K must divide 26; use 2");
    if(bssn::BSSN_RESTORE_SOLVER)
        throw std::runtime_error("Do not restore old/vacuum checkpoints: EKG schema restore checks are not yet wired in this overlay");
    if(bssn::BSSN_GW_EXTRACT_FREQ)
        throw std::runtime_error("Matter-corrected Psi4 extraction is not supplied; set BSSN_GW_EXTRACT_FREQ=0");
    for(unsigned i=0;i<bssn::BSSN_NUM_CONST_VARS_VTU_OUTPUT;++i)
        if(bssn::BSSN_VTU_OUTPUT_CONST_INDICES[i]>3)
            throw std::runtime_error("Only constraint indices 0..3 are computed; Psi4 slots 4,5 are disabled");
    bool has_phi=false,has_pi=false;
    for(unsigned i=0;i<bssn::BSSN_NUM_REFINE_VARS;++i) {
        const auto v=bssn::BSSN_REFINE_VARIABLE_INDICES[i];
        if(v>=26)throw std::runtime_error("Refinement variable index >=26");
        has_phi|=v==24;has_pi|=v==25;
    }
    if(bssn::BSSN_REMESH_TEST_FREQ>0 && (!has_phi || !has_pi))
        throw std::runtime_error("For EKG WAMR include BOTH phi=24 and Pi=25 in REFINE_VARIABLE_INDICES");
    const auto id=bssn::BSSN_ID_TYPE;
    if(id>103)throw std::runtime_error("Unknown EKG initial-data ID; registered new IDs are 100..103");
    if(id==102 && (!EKG_BACKREACTION || !EKG_GAUGE_SYNCHRONOUS || config.mu!=0 ||
                  config.pi0==0 || config.boundary!="analytic"))
        throw std::runtime_error("ID102 FLRW needs coupled+synchronous, mu=0, pi0!=0, analytic boundary");
    if(id==103 && (!EKG_BACKREACTION || EKG_FREEZE_GEOMETRY))
        throw std::runtime_error("ID103 solves the matter Hamiltonian constraint and requires coupled evolution");
    if((id==100 || id==101) && EKG_BACKREACTION && !config.allow_constraint_violating_id)
        throw std::runtime_error("Flat+scalar data are test-field data; use ID103 for a constraint-solved coupled pulse");
    if(config.boundary=="analytic") {
        if(id!=100 && id!=101 && id!=102)throw std::runtime_error("No analytic evolution boundary for this ID");
        if(id!=102 && !EKG_FREEZE_GEOMETRY)throw std::runtime_error("Flat exact boundary requires FREEZE_GEOMETRY=ON");
        if(id==100 && config.mu!=0)throw std::runtime_error("Analytic spherical Gaussian is massless");
        if(EKG_POTENTIAL_AXION && config.mu!=0)throw std::runtime_error("Free plane-wave reference is not an axion exact solution");
    }
    prepare_initial_data(id);
    int rank=0;MPI_Comm_rank(comm,&rank);
    if(!rank) {
        std::cout<<"[EKG] native BSSNCtx/ETS, 26 fields, BACKREACTION="<<EKG_BACKREACTION
                 <<", FREEZE_GEOMETRY="<<EKG_FREEZE_GEOMETRY
                 <<", FUSED_RHS="<<EKG_FUSED_RHS<<"\n";
        if(id==103)std::cout<<"[EKG ID103] ADM mass="<<radial_initial_data().mass
                           <<", Robin residual="<<radial_initial_data().residual<<"\n";
        std::cout<<"[EKG] constraints H, M_i are matter-corrected; Psi4 output is disabled\n";
    }
}
bool initialize_new_id(double x,double y,double z,double* u) {
    const auto id=bssn::BSSN_ID_TYPE;
    if(id<100 || id>103)return false;
    const auto q=initial_state(id,{x,y,z});std::copy(q.begin(),q.end(),u);return true;
}
}
