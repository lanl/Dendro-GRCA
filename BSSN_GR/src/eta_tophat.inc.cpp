// eta_tophat.inc.cpp 
double report = 0.0;
if (dr1 < bhMass1) // near BH1 
  report += bhMass1;
if (dr2 < bhMass2) // near BH2
  report += bhMass2;
if (report == 0.0) // far from both BHs
  report = bhMass1 + bhMass2;
const double eta = 1.0 / report;
