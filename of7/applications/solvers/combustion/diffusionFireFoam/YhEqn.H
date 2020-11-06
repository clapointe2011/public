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
    Qdot = combustion->Qdot();

    volScalarField Yt(0.0*Y[0]);
    
    dEff = turbulence->alphaEff()/(lewisNo + VSMALL);

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
		combustion->R(Yi)
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

    volScalarField& h = thermo.he();

    if (origHEqn)
    {
        fvScalarMatrix hEqn
        (
            fvm::ddt(rho,h) + mvConvection->fvmDiv(phi,h)
          + fvc::ddt(rho,K) + fvc::div(phi,K)
          - fvm::laplacian(turbulence->alphaEff(),h)
          - dpdt
         ==
            combustion->Qdot()
          + radiation->Sh(thermo,h)
          + fvOptions(rho,h)
        );

	hEqn.relax();
        fvOptions.constrain(hEqn);
        hEqn.solve();
    }
    else
    {
        fvScalarMatrix hEqn
        (
            fvm::ddt(rho,h) + mvConvection->fvmDiv(phi,h)
          + fvc::ddt(rho,K) + fvc::div(phi,K)
          + (
                lewisNo != 1
              ? (
		  - fvm::laplacian(turbulence->alphaEff(),h)
                  // based on https://www.cfd-online.com/Forums/openfoam/72127-solving-total-energy-equation.html
		  - (fvc::div(turbulence->muEff()*(fvc::grad(U)&U)) + fvc::div(turbulence->muEff()*(dev2(fvc::grad(U)().T())&U)))
                  // based on https://www.cfd-online.com/Forums/openfoam-programming-development/113264-including-enthalpy-transport.html
                  - fvc::laplacian(h*(dEff-turbulence->alphaEff())*Yt,Yt)
                )
              : (
                  - fvm::laplacian(turbulence->alphaEff(),h)
                  // based on https://www.cfd-online.com/Forums/openfoam/72127-solving-total-energy-equation.html
		  -  (fvc::div(turbulence->muEff()*(fvc::grad(U)&U)) + fvc::div(turbulence->muEff()*(dev2(fvc::grad(U)().T())&U)))
                )
	    )
         ==
            dpdt
          + rho*(U&g)
          + combustion->Qdot()
          + radiation->Sh(thermo,h)
          + fvOptions(rho,h)
        );

        hEqn.relax();
        fvOptions.constrain(hEqn);
        hEqn.solve();
    }

    fvOptions.correct(h);
    thermo.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}
