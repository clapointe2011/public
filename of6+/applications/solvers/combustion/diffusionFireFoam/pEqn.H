rho = thermo.rho();

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

surfaceScalarField phig("phig", -rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
        fvc::flux(rho*HbyA)
      + MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rho, U, phi))
    )
  + phig
);

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phiHbyA, rho, U);
MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf, MRF);

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix p_rghEqn
    (
        fvm::ddt(psi, p_rgh)
      + fvc::ddt(psi, rho)*gh
      + fvc::ddt(psi)*pRef
      + fvc::div(phiHbyA)
      - fvm::laplacian(rhorAUf, p_rgh)
     ==
	    fvOptions(psi, p_rgh, rho.name())
    );

    p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA + p_rghEqn.flux();
        U = HbyA + rAU*fvc::reconstruct((p_rghEqn.flux() + phig)/rhorAUf);
        U.correctBoundaryConditions();
        fvOptions.correct(U);
    }
}

p = p_rgh + rho*gh + pRef;

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

K = 0.5*magSqr(U);

if (mesh.moving())
{
    surfaceVectorField rhoUf = fvc::interpolate(rho*U);
    surfaceVectorField n(mesh.Sf()/mesh.magSf());
    rhoUf += n*(fvc::absolute(phi, rho, U)/mesh.magSf() - (n & rhoUf));
    phi = mesh.Sf() & rhoUf;
}

if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);
    
    if (mesh.moving())
    {
        dpdt -= fvc::div(fvc::meshPhi(rho, U), p);
    }
}
