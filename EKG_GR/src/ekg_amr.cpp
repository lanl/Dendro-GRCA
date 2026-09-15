#include "ekg_native.h"
#include "dataUtils.h"
#include "parameters.h"
#include <algorithm>
#include <array>
#include <utility>
#include <stdexcept>
#include <vector>
namespace ekg {
namespace {
struct RefineView {
    std::array<const double*,26> fields{};
    std::vector<double> phi,pi;
    RefineView(ot::Mesh* mesh,const double** input) {
        std::copy_n(input,26,fields.begin());
        std::size_t n=0;
        // Same block layout as the native unzipped DVector. No cached pointers
        // survive remeshing; every refinement test creates a new read-only view.
        for(const auto& b:mesh->getLocalBlockList())
            n=std::max(n,static_cast<std::size_t>(b.getOffset())+
                static_cast<std::size_t>(b.getAllocationSzX())*b.getAllocationSzY()*b.getAllocationSzZ());
        auto normalize=[&](unsigned f,double scale,std::vector<double>& storage) {
            if(scale==1 || n==0)return;
            storage.resize(n);
            for(std::size_t p=0;p<n;++p)storage[p]=input[f][p]/scale;
            fields[f]=storage.data();
        };
        normalize(24,config.refine_scale_phi,phi);normalize(25,config.refine_scale_pi,pi);
    }
};
}
bool isReMeshWAMR(ot::Mesh* m,const double** u,const unsigned* ids,unsigned n,
                 std::function<double(double,double,double,double*)> tol,double f) {
    RefineView v(m,u);
    return bssn::isReMeshWAMR(m,v.fields.data(),ids,n,std::move(tol),f);
}
bool addRemeshWAMR(ot::Mesh* m,const double** u,const unsigned* ids,unsigned n,
                  std::function<double(double,double,double,double*)> tol,double f,bool relative) {
    RefineView v(m,u);
    if(relative)throw std::runtime_error("Use fixed EKG refinement scales, not native relative_WAMR");
    return bssn::addRemeshWAMR(m,v.fields.data(),ids,n,std::move(tol),f);
}
}
