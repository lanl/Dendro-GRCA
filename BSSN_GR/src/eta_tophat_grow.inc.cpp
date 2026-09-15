// eta_tophat_grow.inc.cpp 
double report = 0.0;
const double t_max = 10.0; // relative time growth stops
const double R1 = .5 * bhMass1 * (1.0 + std::min(t/(t_max*bhMass1), 1.0));
if (dr1 < R1) // near BH1 
  report += bhMass1;
const double R2 = .5 * bhMass2 * (1.0 + std::min(t/(t_max*bhMass2), 1.0));
if (dr2 < R2) // near BH2
  report += bhMass2;
if (report == 0.0) // far from both BHs
  report = bhMass1 + bhMass2;
const double eta = 1.0 / report;
