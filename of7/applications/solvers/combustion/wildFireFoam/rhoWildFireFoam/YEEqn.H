tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);
{   
    combustion->correct();

    volScalarField Yt(0.0*Y[0]);

    dEff = turbulence->alphaEff()/(lewisNo + SMALL);

    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho,Yi)
              + mvConvection->fvmDiv(phi,Yi)
              - fvm::laplacian(dEff,Yi)
             ==
                parcels.SYi(i,Yi)
              + surfaceFilm.Srho(i)
	      + combustion->R(Yi)
              + fvOptions(rho,Yi)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve("Yi");

            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);

    radiation->correct();

    volScalarField& he = thermo.he();
    
    {
        fvScalarMatrix heEqn
        (
            fvm::ddt(rho,he) + mvConvection->fvmDiv(phi,he)
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
            rho*(U&g)
          + combustion->Qdot()
          + radiation->Sh(thermo,he)
          + parcels.Sh(he)
          + surfaceFilm.Sh()
          + fvOptions(rho,he)
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
