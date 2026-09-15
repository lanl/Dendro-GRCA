// eta_single.inc.cpp
const double inpansion = 4.0;
const double rho1 = std::max(dr1, bhMass1/inpansion);
const double rho2 = std::max(dr2, bhMass2/inpansion);
const double rhoG = std::max(bhMass1, bhMass2);
const double eta = 1.0 / std::min({rho1, rho2, rhoG});
