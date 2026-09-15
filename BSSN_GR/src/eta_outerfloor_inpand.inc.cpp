// eta_outerfloor_inpand.inc.cpp
const double inpansion = 16.0; 
const double rho1 = std::max(dr1, bhMass1/inpansion);
const double rho2 = std::max(dr2, bhMass2/inpansion);
const double total_mass = bhMass1 + bhMass2;
// set up radius to enforce 
const double R = 50.0; // M
const double rho_all = std::max(total_mass, r_coord + total_mass - R);
const double eta = 1.0 / std::min({rho1, rho2, rho_all});
