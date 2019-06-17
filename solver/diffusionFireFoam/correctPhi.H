{
    const fvMesh& mesh = U.mesh();
    const Time& runTime = mesh.time();

    correctUphiBCs(rho, U, phi);
    correctUphiBCs(rho.oldTime(), U.oldTime(), phi.oldTime());

    // Initialize BCs list for pcorr to zero-gradient
    wordList pcorrTypes
    (
        p_rgh.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    // Set BCs of pcorr to fixed-value for patches at which p is fixed
    forAll(p_rgh.boundaryField(), patchi)
    {
        if (p_rgh.boundaryField()[patchi].fixesValue())
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
    
    volScalarField pcorr_0
    (
        IOobject
        (
            "pcorr_0",
            mesh
        ),
        mesh,
        dimensionedScalar("pcorr_0", p.dimensions(), 0.0),
        pcorrTypes
    );

    mesh.setFluxRequired(pcorr.name());
    mesh.setFluxRequired(pcorr_0.name());
    
    dimensionedScalar rAUf("rAUf", dimTime, 1.0);

    //CL flux corrections
    for (int i = 0; i<nCorrPhi; i++)
    {
        // Solve for pcorr such that the divergence of the corrected flux
        // matches the divRhoU provided (from previous iteration, time-step...)
        fvScalarMatrix pcorrEqn
        (
            fvc::div(phi)
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
        
        // Solve for pcorr such that the divergence of the corrected flux
        // matches the divRhoU provided (from previous iteration, time-step...)
        fvScalarMatrix pcorrEqn_0
        (
            fvc::div(phi.oldTime())
          - fvm::laplacian(rAUf, pcorr_0)
         ==
            divrhoU_0()
        );

        pcorrEqn_0.solve
        (
            mesh.solver(pcorr.select(pimple.finalNonOrthogonalIter()))
        );

        if (pimple.finalNonOrthogonalIter())
        {
            phi.oldTime() += pcorrEqn_0.flux();
        }
    }
    
    if (ddtSchemeName == "backward")
    {
        correctUphiBCs(rho.oldTime().oldTime(), U.oldTime().oldTime(), phi.oldTime().oldTime());
        
        volScalarField pcorr_0_0
        (
            IOobject
            (
                "pcorr_0_0",
                mesh
            ),
            mesh,
            dimensionedScalar("pcorr_0_0", p.dimensions(), 0.0),
            pcorrTypes
        );
        
        mesh.setFluxRequired(pcorr_0_0.name());
        
        //CL flux corrections
        for (int i = 0; i<nCorrPhi; i++)
        {   
            // Solve for pcorr such that the divergence of the corrected flux
            // matches the divRhoU provided (from previous iteration, time-step...)
            fvScalarMatrix pcorrEqn_0_0
            (
                fvc::div(phi.oldTime().oldTime())
              - fvm::laplacian(rAUf, pcorr_0_0)
             ==
                divrhoU_0_0()
            );

            pcorrEqn_0_0.solve
            (
                mesh.solver(pcorr.select(pimple.finalNonOrthogonalIter()))
            );

            if (pimple.finalNonOrthogonalIter())
            {
                phi.oldTime().oldTime() += pcorrEqn_0_0.flux();
            }
        }
    }
}
