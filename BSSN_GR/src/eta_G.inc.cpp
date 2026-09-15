// eta_G.inc.cpp
// extreme mass ratio shift damper eta_G from
// [Lousto & Healy '23](https://arxiv.org/abs/2203.08831)

// define useful constants:
// distance from each BH to grid center
const double r1 = sqrt(bh1x*bh1x + bh1y*bh1y + bh1z*bh1z);
const double r2 = sqrt(bh2x*bh2x + bh2y*bh2y + bh2z*bh2z);
// exponential softening scales
const double s1 = 2.0 * bhMass1;
const double s2 = 2.0 * bhMass2;
// Hardcoding A = B = C = 1, n = 2 as in LH23 formulation
const double eta = 1.0 / (bhMass1 + bhMass2)
                 + (1.0 / bhMass1) * pow(r1 * r1 / (r1 * r1 + s2 * s2), 2) * exp(-(dr1 * dr1) / (s1 * s1))
                 + (1.0 / bhMass2) * pow(r2 * r2 / (r2 * r2 + s1 * s1), 2) * exp(-(dr2 * dr2) / (s2 * s2));

