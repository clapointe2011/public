Info << "PISO iteration: " << pimple.corrPISO() << endl;

if (pimple.corrPISO() > 1)
{
    thermo.correct();
}

//CL bound density
thermo.rho() = max(thermo.rho(), rhoMin);
thermo.rho() = min(thermo.rho(), rhoMax);
thermo.rho().relax();

thermo.rho() -= psi*p;
rhoHat = thermo.rho();

Info << "max/min of rhoHat :" << max(rhoHat).value() << "/" << min(rhoHat).value() << endl;

volVectorField HbyA ("HbyA", U);
volScalarField rAU ("rAU", 1.0 / UEqn.A());
HbyA = UEqn.H() * rAU;

psi_pos = fvc::interpolate(psi, pos, "reconstruct(psi)");
psi_neg = fvc::interpolate(psi, neg, "reconstruct(psi)");

psiU_pos= fvc::interpolate(psi*HbyA, pos, "reconstruct(U)");
psiU_neg= fvc::interpolate(psi*HbyA, neg, "reconstruct(U)");

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
phid_pos = aphiv_pos * psi_pos;
phid_neg = aphiv_neg * psi_neg;

rhoHat_pos=fvc::interpolate(rhoHat, pos, "reconstruct(rhoHat)");
rhoHat_neg=fvc::interpolate(rhoHat, neg, "reconstruct(rhoHat)");

rhoHatPhi_pos = aphiv_pos * rhoHat_pos;
rhoHatPhi_neg = aphiv_neg * rhoHat_neg;

surfaceScalarField Dp_pos
(
    "Dp_pos",
    fvc::interpolate(rho*rAU, pos, "reconstruct(Dp)")
);

surfaceScalarField Dp_neg
(
    "Dp_neg",
    fvc::interpolate(rho*rAU, neg, "reconstruct(Dp)")
);

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn_pos
    (
	    fvc::div(rhoHatPhi_pos) + fvm::div(phid_pos,p) - fvm::laplacian(Dp_pos*a_pos,p)
    );
    
    fvScalarMatrix pEqn_neg
    (
	    fvc::div(rhoHatPhi_neg) + fvm::div(phid_neg,p) - fvm::laplacian(Dp_neg*a_neg,p)
    );
    
    fvScalarMatrix pEqn
    (
          fvm::ddt(psi,p)
        + fvc::ddt(rhoHat)
        + pEqn_pos
        + pEqn_neg
        ==
          fvOptions(psi,p,rho.name())
    );
    
    fvOptions.constrain(pEqn);
    
    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));
    
    if (pimple.finalNonOrthogonalIter())
    {
        phiPos = pEqn_pos.flux() + rhoHatPhi_pos;
        phiNeg = pEqn_neg.flux() + rhoHatPhi_neg;
        phi = phiPos + phiNeg;
    }
}

p_pos = fvc::interpolate(p, pos, "reconstruct(p)");
p_neg = fvc::interpolate(p, neg, "reconstruct(p)");

gradp = fvc::div((a_pos*p_pos + a_neg*p_neg)*mesh.Sf());
#include "limitPressureGradient.H"
U = HbyA - rAU * gradp;
U.correctBoundaryConditions();
fvOptions.correct(U);
K = 0.5*magSqr(U);

Info << "max(U): " << max(U).value() << endl;
Info << "max/min of p: " << max(p).value() << "/" << min(p).value() << endl;

#include "rhoEqn.H"
thermo.rho() += psi*p;
#include "compressibleContinuityErrs.H"

rho = thermo.rho();
rhof = p_pos * psi_pos * a_pos + rhoHat_pos*a_pos + p_neg * psi_neg * a_neg + rhoHat_neg*a_neg;

if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);
}
