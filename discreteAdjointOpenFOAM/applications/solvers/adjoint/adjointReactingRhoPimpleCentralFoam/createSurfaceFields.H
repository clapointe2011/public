surfaceScalarField psi_pos
(
    "psi_pos",
    fvc::interpolate(psi, pos, "reconstruct(psi)")
);

surfaceScalarField psi_neg
(
    "psi_neg",
    fvc::interpolate(psi, neg, "reconstruct(psi)")
);

surfaceVectorField psiU_pos
(
    "psiU_pos",
    fvc::interpolate(psi*U, pos, "reconstruct(U)")
);

surfaceVectorField psiU_neg
(
    "psiU_neg",
    fvc::interpolate(psi*U, neg, "reconstruct(U)")
);

surfaceScalarField p_pos
(
    "p_pos",
    fvc::interpolate(p, pos, "reconstruct(p)")
);

surfaceScalarField p_neg
(
    "p_neg",
    fvc::interpolate(p, neg, "reconstruct(p)")
);
    
surfaceScalarField phiv_pos
(
    "phiv_pos",
    (psiU_pos / psi_pos)  & mesh.Sf()
);

surfaceScalarField phiv_neg
(
    "phiv_neg",
    (psiU_neg / psi_neg) & mesh.Sf()
);

volScalarField c
(
    //"cSound",
    //sqrt(thermo.Cp()/thermo.Cv() / psi)

    IOobject
    (
        "cSound",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    sqrt(thermo.Cp()/thermo.Cv() / psi)
);
    
surfaceScalarField cSf_pos
(
    "cSf_pos",
    fvc::interpolate(c, pos, "reconstruct(psi)")*mesh.magSf()
);

surfaceScalarField cSf_neg
(
    "cSf_neg",
    fvc::interpolate(c, neg, "reconstruct(psi)")*mesh.magSf()
);

surfaceScalarField ap
(
    "ap",
    max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero)
);
surfaceScalarField am
(
    "am",
    min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero)
);

surfaceScalarField a_pos
(
    "a_pos",
    ap/(ap - am)
);

surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));

surfaceScalarField aSf
(
    "aSf",
    am*a_pos
);

surfaceScalarField a_neg
(
    "a_neg",
    1.0 - a_pos
);

phiv_pos *= a_pos;
phiv_neg *= a_neg;

surfaceScalarField aphiv_pos("aphiv_pos", phiv_pos - aSf);
surfaceScalarField aphiv_neg("aphiv_neg", phiv_neg + aSf);

surfaceScalarField phid_pos
(
    "phid_pos",
    aphiv_pos * psi_pos
);

surfaceScalarField phid_neg
(
    "phid_neg",
    aphiv_neg * psi_neg
);

surfaceScalarField phiPos
(
    "phiPos",
    phi
);

surfaceScalarField phiNeg
(
    "phiNeg",
    phi * 0.0
);

surfaceScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimless
);

surfaceScalarField rhoHat_pos
(
    "rhoHat_pos",
    fvc::interpolate(rhoHat, pos, "reconstruct(rhoHat)")
);

surfaceScalarField rhoHat_neg
(
    "rhoHat_neg",
    fvc::interpolate(rhoHat, neg, "reconstruct(rhoHat)")
);

surfaceScalarField rhoHatPhi_pos
(
    "rhoHatPhi_pos",
    phiPos*0.0
);

surfaceScalarField rhoHatPhi_neg
(
     "rhoHatPhi_neg",
     phiNeg*0.0
);
    
surfaceScalarField rhof
(
    "rhof",
    p_pos*psi_pos*a_pos + rhoHat_pos*a_pos + p_neg*psi_neg*a_neg + rhoHat_neg*a_neg
);

forAll(kappa, iFace)
{
    kappa[iFace] = 0.0;
}
    
forAll(kappa.boundaryField(), iPatch)
{
    kappa.boundaryField()[iPatch] = 0.0;
}
