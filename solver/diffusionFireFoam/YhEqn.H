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
    radiation->correct();
    combustion->correct();
    Qdot = combustion->Qdot();
    volScalarField Yt(0.0*Y[0]);
    
    dEff = turbulence->alphaEff()/(lewisNo + SMALL);

    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(dEff, Yi)
             ==
		        combustion->R(Yi)
              + fvOptions(rho, Yi)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve(mesh.solver("Yi"));

            fvOptions.correct(Yi);

            Yi.max(0.0);
            Yt += Yi;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);

    volScalarField& h = thermo.he();

    fvScalarMatrix hEqn
    (
        fvm::ddt(rho, h) + mvConvection->fvmDiv(phi, h)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      - fvm::laplacian(turbulence->alphaEff(), h)
      - dpdt
     ==
        Qdot
      + radiation->Sh(thermo, h)
      + fvOptions(rho, h)
    );

    hEqn.relax();

    fvOptions.constrain(hEqn);

    hEqn.solve();

    fvOptions.correct(h);

    thermo.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}
