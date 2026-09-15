#pragma once
#include "ekg_config.h"
#include "ekg_build_config.h"
#include "ekg_generated_config.h"
#include <functional>
#include <mpi.h>
#include <toml.hpp>
namespace ot { class Mesh; }
namespace ekg {
// Zero disables optional native events without modulo-by-zero.
inline bool event_due(unsigned long long step,unsigned long long frequency) {
    return frequency != 0 && step % frequency == 0;
}
void read_config(const toml::value& root);
toml::value config_as_toml();
void validate_runtime(unsigned ts_mode,MPI_Comm comm);
bool initialize_new_id(double x,double y,double z,double* state);
void rhs_block(double** out,const double** in,unsigned offset,const double* pmin,
               const double* pmax,const unsigned* size,unsigned bflag,double stage_time);
void constraints_block(double** out,const double** in,unsigned offset,
                       const double* pmin,const double* pmax,const unsigned* size,unsigned bflag);
bool isReMeshWAMR(ot::Mesh* mesh,const double** values,const unsigned* ids,unsigned n,
                 std::function<double(double,double,double,double*)> tolerance,double coarsen);
bool addRemeshWAMR(ot::Mesh* mesh,const double** values,const unsigned* ids,unsigned n,
                  std::function<double(double,double,double,double*)> tolerance,double coarsen,
                  bool relative=false);
}
