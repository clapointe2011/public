MRF.correctBoundaryVelocity(U);

fvVectorMatrix UEqn
(
      fvm::ddt(rho,U)
    + fvm::div(phiPos,U)
    + fvm::div(phiNeg,U)
    + rho*fvm::Sp(PRweight*alpha, U)
    + MRF.DDt(rho,U)
    + turbulence->divDevRhoReff(U)
    ==
      fvOptions(rho, U)
);

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve
    (
	    UEqn == -gradp
    );
    
    fvOptions.correct(U);
}
