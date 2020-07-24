Info<< "\nStarting Final Loop\n" << endl;

//CL final passive loop

for(int fIter = 1; fIter <= finalCutoff; fIter++){

	Info << "\npassive solution time : " << fIter << nl << endl;

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
		    turbulence->correct();
	    }
	}
}
