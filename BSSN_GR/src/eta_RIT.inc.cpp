// RIT constant eta formulation
const double w        = r_coord / bssn::RIT_ETA_WIDTH;
const double arg      = -w * w * w * w;
const double eta = (bssn::RIT_ETA_CENTRAL - bssn::RIT_ETA_OUTER) * exp(arg) + bssn::RIT_ETA_OUTER;

