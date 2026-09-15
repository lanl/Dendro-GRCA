// eta_global_ramp.inc.cpp
constexpr double eta_init   = 4.00;   // strong early damping of gauge waves
constexpr double eta_final  = 2.00;   // your proven good late-time value
constexpr double tau        = 10.0;   // decay timescale in M; tune 15–40
const double eta = eta_final + (eta_init - eta_final) * std::exp(-t / tau);
