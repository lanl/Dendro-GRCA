// eta_causal_grow.inc.cpp

// minimum distance to a given BH's position;
// technically should maybe be start position.
const double r_min = std::min(dr1, dr2);

// eta damps shift vector oscillations; we want: 
//  - small value in timelike regions, for freer motion
//  - large value in spacelike regions, to damp the initial gauge wave
constexpr double ETA_TIMELIKE  = 0.25;  // free/fast motion regime
constexpr double ETA_SPACELIKE = 2.50;  // strong damping on gauge waves

// Select eta based on whether the point is causally connected
const double eta = (r_min < t) ? ETA_TIMELIKE : ETA_SPACELIKE;
