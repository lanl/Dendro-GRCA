/**
 * @file bssn_gw_adapter.h
 * @brief Maps BSSN's GW parameters onto dendrolib's shared GW extractor.
 *
 * `dendro_gr::extractFarFieldPsi4` (dendrolib GR/include/gw_extract.h) is the
 * theory-agnostic version of what BSSN_GR/include/gwExtract.h has always done
 * locally: it takes a GWExtractionConfig instead of reading BSSN globals, so
 * the same code serves BSSN, CCZ4, Z4c and EMDA. This header is the only place
 * that knows how BSSN's parameters map onto that config.
 *
 * Selected by BSSN_GW_USE_DENDROLIB (default ON). With the option OFF the
 * BSSN-local GW::extractFarFieldPsi4 is used instead and this header is not
 * compiled.
 */

#ifndef BSSN_GW_ADAPTER_H
#define BSSN_GW_ADAPTER_H

#include "grDef.h"
#include "gw_extract.h"
#include "parameters.h"

namespace bssn {

/**
 * @brief Build the extraction config from the currently loaded BSSN parameters.
 *
 * Cheap enough to call per extraction (a few small vector copies), and doing so
 * keeps it correct if a parameter is ever changed mid-run.
 */
inline dendro_gr::GWExtractionConfig makeGWExtractionConfig() {
    dendro_gr::GWExtractionConfig cfg;

    cfg.psi4_real_idx = bssn::VAR_CONSTRAINT::C_PSI4_REAL;
    cfg.psi4_imag_idx = bssn::VAR_CONSTRAINT::C_PSI4_IMG;

    // radii stay double -- par files parse them with as_floating()
    cfg.radii.assign(GW::BSSN_GW_RADAII,
                     GW::BSSN_GW_RADAII + GW::BSSN_GW_NUM_RADAII);
    cfg.l_modes.assign(GW::BSSN_GW_L_MODES,
                       GW::BSSN_GW_L_MODES + GW::BSSN_GW_NUM_LMODES);

    cfg.file_prefix = bssn::BSSN_PROFILE_FILE_PREFIX;

    for (unsigned int d = 0; d < 3; d++) {
        cfg.compd_min[d]  = bssn::BSSN_COMPD_MIN[d];
        cfg.compd_max[d]  = bssn::BSSN_COMPD_MAX[d];
        cfg.octree_min[d] = bssn::BSSN_OCTREE_MIN[d];
        cfg.octree_max[d] = bssn::BSSN_OCTREE_MAX[d];
    }

    cfg.output_precision = GW::BSSN_GW_OUTPUT_PRECISION;

    return cfg;
}

}  // namespace bssn

#endif  // BSSN_GW_ADAPTER_H
