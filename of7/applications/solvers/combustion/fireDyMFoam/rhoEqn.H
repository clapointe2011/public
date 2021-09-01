{
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(rho)
      + fvc::div(phi)
      ==
        parcels.Srho(rho)
      + surfaceFilm.Srho()
	  + fvOptions(rho)
    );

    rhoEqn.solve();

    fvOptions.correct(rho);
}
