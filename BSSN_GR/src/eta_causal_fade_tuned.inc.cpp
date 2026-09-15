// eta_causal_fade_tuned.inc.cpp
constexpr double ETA_SPACELIKE = 3.50;   // stronger early damping
constexpr double ETA_TIMELIKE  = 1.75;   // late-time free motion
constexpr double T_FADE        = 50.0;   // slower = safer (was 20)
const double t_star = std::clamp(t / T_FADE, 0.0, 1.0);  // global time, not (t-r_min)
const double eta = std::pow(ETA_SPACELIKE, 1.0 - t_star) * std::pow(ETA_TIMELIKE, t_star);  // geometric is smoothest
