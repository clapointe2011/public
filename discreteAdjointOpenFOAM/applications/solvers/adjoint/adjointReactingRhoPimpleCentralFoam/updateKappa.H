aphiv_pos = phiPos / (p_pos * psi_pos);
aphiv_neg = phiNeg / (p_neg * psi_neg);
amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));

surfaceScalarField amaxSfbyDelta
(
    mesh.surfaceInterpolation::deltaCoeffs()*amaxSf
);

surfaceScalarField Maf
(
    "Maf",
    mag(phi) / (psi_pos*p_pos*a_pos + psi_neg*p_neg*a_neg)
    / (cSf_pos*a_pos + cSf_neg*a_neg)
);

Info << "max/min Maf: " << max(Maf).value() << "/" << min(Maf).value() << endl;

dimensionedScalar cDeltaT = runTime.deltaT();

if (LTS)
{
    cDeltaT.value() = 1.0 / gMax(trDeltaT().internalField());
}

kappa = 
    max
    (
        min
        (
            Maf / (amaxSfbyDelta/mesh.magSf() * cDeltaT),
            scalar(1.0)
        ),
        scalar(0.0)
    );

forAll(kappa.boundaryField(), iPatch)
{
    fvsPatchField<scalar>& kappapf = kappa.boundaryField()[iPatch];
    if (isA<coupledFvsPatchField<scalar> > (kappapf))
    {
	    forAll (kappapf, iFace)
	    {
	        kappapf[iFace] = 0.0;
	    }
    }
}

Info << "max / min kappa: " << max(kappa).value() << "/" << min(kappa).value() << endl;
phiPos = phiPos + (1.0 - kappa) * phiNeg;
phiNeg = kappa * phiNeg;


