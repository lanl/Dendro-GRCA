// eta_causal_fade.inc.cpp

// minimum distance to a given BH's position;
// technically should maybe be start position.
const double r_min = std::min(dr1, dr2);

// eta damps shift vector oscillations; we want: 
// - small value in timelike regions, for freer motion
// - large value in spacelike regions, to damp the initial gauge wave
constexpr double ETA_TIMELIKE  = 1.00; // free/fast motion regime
constexpr double ETA_SPACELIKE = 2.00; // strong damping on gauge waves

// time to fade from one to the other
// should be ok for lower mass ratios
constexpr double T_FADE = 20.00; 

// relative time variable, from 0 to 1
const double t_star = std::clamp((t - r_min) / T_FADE, 0.0, 1.0);

// Select eta based on whether the point is causally connected
// arithmetic slide: 
// const double eta = (1.0 - t_star)*ETA_SPACELIKE + t_star*ETA_TIMELIKE;
// geometric slide:
const double eta = std::pow(ETA_SPACELIKE, 1.0 - t_star) * std::pow(ETA_TIMELIKE, t_star);

