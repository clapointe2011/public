fvScalarMatrix hEqn
(
      fvm::ddt(rho, h) + fvm::div(phiPos, h) + fvm::div(phiNeg, h)
    + fvc::ddt(rho, K) + fvc::div(phiPos, K) + fvc::div(phiNeg, K)
    - fvm::laplacian(turbulence->alphaEff(), h)
    //- fvc::div( ((-turbulence->devRhoReff()) & U) )
    - dpdt
    ==
      reaction->Sh()
    + fvOptions(rho,h)
);

hEqn.relax();

fvOptions.constrain(hEqn);

hEqn.solve();

fvOptions.correct(h);

thermo.correct();

Info << "max/min T:" << max(T).value() << "/" << min(T).value() << endl;
