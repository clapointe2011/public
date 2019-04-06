{
    const fvMesh& mesh = U.mesh();
    const Time& runTime = mesh.time();

    correctUphiBCs(rho, U, phi);

    // Initialize BCs list for pcorr to zero-gradient
    wordList pcorrTypes
    (
        p.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    // Set BCs of pcorr to fixed-value for patches at which p is fixed
    forAll(p.boundaryField(), patchi)
    {
        if (p.boundaryField()[patchi].fixesValue())
        {
            pcorrTypes[patchi] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("pcorr", p.dimensions(), 0.0),
        pcorrTypes
    );
    
    dimensionedScalar rAUf("rAUf", dimTime, 1.0);
    //surfaceScalarField rhorAUf("rhorAUf",fvc::interpolate(rho*rAU()));

    mesh.setFluxRequired(pcorr.name());

    while (pimple.correct())
    {
        // Solve for pcorr such that the divergence of the corrected flux
        // matches the divRhoU provided (from previous iteration, time-step...)
        fvScalarMatrix pcorrEqn
        (
            fvm::ddt(psi, pcorr)
          + fvc::div(phi)
          - fvm::laplacian(rAUf, pcorr)
         ==
            divrhoU()
        );

        pcorrEqn.solve
        (
            mesh.solver(pcorr.select(pimple.finalNonOrthogonalIter()))
        );

        if (pimple.finalNonOrthogonalIter())
        {
            phi += pcorrEqn.flux();
        }
    }
}
