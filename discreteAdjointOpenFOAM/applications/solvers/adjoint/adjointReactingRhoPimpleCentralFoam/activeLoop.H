Info<< "\nactive iteration : " << aIter << nl << endl;

if (LTS)
{
    #include "setRDeltaT.H"
}
else
{
    #include "acousticCourantNo.H"
    #include "compressibleCourantNo.H"
    #include "readTimeControls.H"
    #include "setDeltaT.H"
}
        
runTime++;
        
/*rho.oldTime();
rhoHat.oldTime();
p.oldTime();
U.oldTime();
h.oldTime();
Ek.oldTime();
psi.oldTime();

forAll(Y, i){
    Y[i].oldTime();
}*/

#include "rhoEqn.H"

//PIMPLE loop
while (pimple.loop())
{
    #include "UEqn.H"
    #include "YEqn.H"
    #include "hEqn.H"

    //Pressure Corrector
    while (pimple.correct())
    {
	    #include "pEqn.H"
    }

    #include "updateKappa.H"

    if (pimple.turbCorr())
    {
	    //CL allow frozen turbulence
	    if (frozenTurbulence){
	        dco::ga1s<double>::global_tape->switch_to_passive();
	    }

	    turbulence->correct();

	    if (frozenTurbulence){
	        dco::ga1s<double>::global_tape->switch_to_active();
	    }
    }
}
