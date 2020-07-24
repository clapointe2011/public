/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    pisoCentralFoam

Description
    Pressure-based semi implicit compressible flow solver based on central-upwind schemes of
    Kurganov and Tadmor

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rhoThermo.H"
#include "wallFvPatch.H"
#include "pimpleControl.H"
#include "rhoCombustionModel.H"
#include "turbulentFluidThermoModel.H"

#include "fixedFluxPressureFvPatchScalarField.H"
#include "coupledFvsPatchFields.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "cellQuality.H"
#include "fvIOoptionList.H"

#include "costFunctionLibrary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    
    #include "createTime.H"
    #include "createMesh.H"
    
    pimpleControl pimple(mesh);
    
    #include "createTimeControls.H"
    #include "createRDeltaT.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "createMRF.H"
    #include "initContinuityErrs.H"
    #include "readCourantType.H"
    
    #include "createSurfaceFields.H"
    #include "markBadQualityCells.H"
    
    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    Info<< "\nStarting Optimization\n" << endl;
    
    #include "initKappaField.H"
    
    dco::ga1s<double>::global_tape = dco::ga1s<double>::tape_t::create();
    
    //CL passive loop
	if (passiveInitialization){
	    dco::ga1s<double>::global_tape->switch_to_passive();
	    #include "passiveInitialization.H"
	}

    for (int nIter = 1; nIter <= optEndIter; nIter++){

	    Info << "optimization iteration : " << nIter << nl << endl;
	    
	    //CL clear tape
	    dco::ga1s<double>::global_tape->reset(); //CL

	    //CL switch tape to active
	    dco::ga1s<double>::global_tape->switch_to_active();

	    //CL register alpha
	    forAll(alpha,i){
	        dco::ga1s<double>::global_tape->register_variable(alpha[i]);
	    }

	    //CL active loop
	    for(int aIter = 1; aIter <= tapedIter; aIter++){
	        #include "activeLoop.H"
	    }

	    //CL calculate updated cost function
	    Info << "\ncalculating updated cost function" << endl; //CL
	    scalar J = CostFunction(mesh).eval(); //CL
	    Info << "new cost function value : " << J << nl << endl; //CL

	    // Set adjoint of J
	    if (Pstream::master()){
	        dco::derivative(J)=1;
	    }

	    //CL interpret adjoints
	    Info<< "interpreting tape" << endl;
	    dco::ga1s<double>::global_tape->interpret_adjoint();	

	    //CL calculate sensitivities
	    #include "calcSens.H"

	    //CL switch to passive to be safe
	    dco::ga1s<double>::global_tape->switch_to_passive();

	    //CL update alpha
	    alpha += mesh.fieldRelaxationFactor("alpha")*(min(max(alpha + lambda*(-sens), zeroAlpha), alphaMax) - alpha);

	    //CL process alpha field
	    #include "procAlpha.H"

	    //CL write fields for iteration nIter
	    if (writeOptFields){
	        if (nIter % writeOptFreq == 0){
		        runTime.setTime(nIter,(label)(nIter)); //CL
		        runTime.writeNow(); //CL
	        }
	    }

	    //CL write alpha to 0 directory
	    if (cumulativeSolution){
	        runTime.setTime(0,(label)0); //CL
	        runTime.writeNow(); //CL
	    }else{
	        runTime.setTime(0,(label)0); //CL
	        alpha.write(); //CL
	    }

	    //CL check for convergence
	    #include "solutionControl.H" //CL
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
