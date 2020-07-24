psi_pos = fvc::interpolate(psi, pos, "reconstruct(psi)");
psi_neg = fvc::interpolate(psi, neg, "reconstruct(psi)");

p_pos = fvc::interpolate(p, pos, "reconstruct(p)");
p_neg = fvc::interpolate(p, neg, "reconstruct(p)");

psiU_pos= fvc::interpolate(psi*U, pos, "reconstruct(U)");
psiU_neg= fvc::interpolate(psi*U, neg, "reconstruct(U)");

phiv_pos= (psiU_pos / psi_pos) & mesh.Sf();
phiv_neg= (psiU_neg / psi_neg) & mesh.Sf();

MRF.makeRelative(phiv_pos);
MRF.makeRelative(phiv_neg);

c = sqrt(thermo.Cp()/thermo.Cv() / psi);

cSf_pos = fvc::interpolate(c, pos, "reconstruct(psi)")*mesh.magSf();
cSf_neg = fvc::interpolate(c, neg, "reconstruct(psi)")*mesh.magSf();

ap = max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero);
am = min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero);

a_pos = ap/(ap - am);
aSf = am*a_pos;
a_neg = 1.0 - a_pos;

phiv_pos *= a_pos;
phiv_neg *= a_neg;
aphiv_pos = phiv_pos - aSf;
aphiv_neg = phiv_neg + aSf;

phiPos = aphiv_pos * p_pos * psi_pos;
phiNeg = aphiv_neg * p_neg * psi_neg;
phi = phiPos + phiNeg;

#include "updateKappa.H"

