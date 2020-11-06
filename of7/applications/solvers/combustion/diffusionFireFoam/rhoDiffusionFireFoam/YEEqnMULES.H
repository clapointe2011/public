//CL MULES implementation based on similar functionality in
//reactingPimpleCentralFoam from Kraposhin et. al.
//and foam solver interMixingFoam

{
    combustion->correct();

    volScalarField Yt(0.0*Y[0]);
    
    dEff = turbulence->alphaEff()/(lewisNo + VSMALL);

    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix SuSp
            (
                combustion->R(Yi)
              + fvOptions(rho,Yi)
            );
            volScalarField Su = SuSp & Yi;

	    const word YiName(Yi.name());
            Yi.rename("Yi");

            //Bounded (implicit) Flux
            fvScalarMatrix divPhiYi
            (
                fv::gaussConvectionScheme<scalar>
                (
                    mesh,
                    phi,
                    upwind<scalar>(mesh, phi)
                ).fvmDiv(phi, Yi)
            );

            //Correction Flux
            surfaceScalarField phiYiCorr
            (
                fvc::flux
                (
                    phi,
                    Yi,
                    "div(phi,Yi)"
                )
              - divPhiYi.flux()
            );

            MULES::limitCorr
            (
	        1/mesh.time().deltaTValue(),
                rho,
                Yi,
                divPhiYi.flux(),
                phiYiCorr,
                zeroField(),
                Su,
                oneField(),
                zeroField()
            );
            Yi.rename(YiName);

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho,Yi)
              + divPhiYi + fvc::div(phiYiCorr)
              - fvm::laplacian(dEff, Yi)
	      ==
                SuSp
            );

            YiEqn.relax();
            fvOptions.constrain(YiEqn);
            YiEqn.solve("Yi");
            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;
        }
    }

    Y[inertIndex] = scalar(1.0) - Yt;
    Y[inertIndex].max(0.0);

    radiation->correct();

    volScalarField& he = thermo.he();

    //Bounded (implict) Flux
    fvScalarMatrix divPhiHE
    (
        fv::gaussConvectionScheme<scalar>
        (
            mesh,
            phi,
            upwind<scalar>(mesh, phi)
        ).fvmDiv(phi, he)
    );

    //Correction Flux
    surfaceScalarField phiHECorr
    (
        fvc::flux
        (
            phi,
            he,
            "div(phi,he)"
        )
      - divPhiHE.flux()
    );

    fvScalarMatrix SuSp
    (
        combustion->Qdot()
      + rho*(U&g)
      + radiation->Sh(thermo,he)
      + fvOptions(rho,he)
    );
    volScalarField Su = SuSp & he;

    scalar Tmin = pimple.dict().lookupOrDefault<scalar>("Tmin",300);
    scalar Tmax = pimple.dict().lookupOrDefault<scalar>("Tmax",3000);
    scalarField heSmall = thermo.Cp()*Tmin;
    scalarField heLarge = thermo.Cp()*Tmax;

    MULES::limitCorr
    (
        1/mesh.time().deltaTValue(),
        rho,
        he,
        divPhiHE.flux(),
        phiHECorr,
        zeroField(),
        Su,
        heLarge,
        heSmall
    );

    {
        fvScalarMatrix heEqn
        (
            fvm::ddt(rho,he) + divPhiHE
          + fvc::div(phiHECorr)
          + fvc::ddt(rho,K) + fvc::div(phi,K)
          + (
                he.name() == "e"
              ? fvc::div
                (
                    fvc::absolute(phi/fvc::interpolate(rho), U),
                    p,
                    "div(phiv,p)"
                )
                : -dpdt
            )
          - fvm::laplacian(turbulence->alphaEff(),he)
          ==
            SuSp
        );

        heEqn.relax();
        fvOptions.constrain(heEqn);
        heEqn.solve();
        fvOptions.correct(he);
    }
    
    thermo.correct();

    //Update stored Qdot after Yi/he solve
    Qdot = combustion->Qdot();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}
