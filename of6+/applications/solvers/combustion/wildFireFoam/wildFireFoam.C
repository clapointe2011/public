/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    fireFoam

Description
    Transient solver for fires and turbulent diffusion flames with reacting
    particle clouds, surface film and pyrolysis modelling.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "basicReactingMultiphaseCloud.H"
#include "surfaceFilmModel.H"
#include "pyrolysisModelCollection.H"
#include "radiationModel.H"
#include "SLGThermo.H"
#include "solidChemistryModel.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "dynamicFvMesh.H"
#include "motionSolver.H"
#include "myDynamicRefineFvMesh.H"
#include "fvcSmooth.H"
#include "cellSet.H"
#include "CorrectPhi.H"
#include "localEulerDdtScheme.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createCustomControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    
    autoPtr<volVectorField> rhoU;
    {
        rhoU = new volVectorField("rhoU", rho*U);
    }
    
    autoPtr<volVectorField> rhoU_0;
    {
        rhoU_0 = new volVectorField("rhoU_0", rho*U);
    }

    autoPtr<volVectorField> rhoU_0_0;
    {
        rhoU_0_0 = new volVectorField("rhoU_0_0", rhoU);
    }
    
    word ddtSchemeName(U.mesh().ddtScheme(U.name()));

    // if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }
    
    #include "readPyrolysisTimeControls.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    
    startRunIndex = runTime.timeIndex();

    while (runTime.run())
    {
	    #include "readCustomControls.H"
	    
        autoPtr<volScalarField> divrhoU;
        if (correctPhi)
        {
            divrhoU = new volScalarField
            (
                "divrhoU",
                fvc::div(fvc::absolute(phi, rho, U))
            );
        }
        
        autoPtr<volScalarField> divrhoU_0;
        if (correctPhi)
        {
            divrhoU_0 = new volScalarField
            (
                "divrhoU_0",
                fvc::div(fvc::absolute(phi.oldTime(), rho.oldTime(), U.oldTime()))
            );
        }
        
        autoPtr<volScalarField> divrhoU_0_0;
        if (correctPhi && ddtSchemeName != "Euler")
        {
            divrhoU_0_0 = new volScalarField
            (
                "divrhoU_0_0",
                fvc::div(fvc::absolute(phi.oldTime().oldTime(), rho.oldTime().oldTime(), U.oldTime().oldTime()))
            );
        }
        
        if (solvePrimaryRegion && (runTime.timeIndex() - startRunIndex) >= startDyMIndex)
        {
            #include "updateMesh.H"
        }

        /*if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else*/
        {
            #include "calcCourantNo.H"
            #include "solidRegionDiffusionNo.H"
            #include "setMultiRegionDeltaT.H"
            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        parcels.evolve();

        surfaceFilm.evolve();

        if (solvePyrolysisRegion)
        {
            pyrolysis.evolve();
        }

        if (solvePrimaryRegion)
        {
        
            #include "rhoEqn.H"

            // --- PIMPLE loop
            while (pimple.loop())
            {
                #include "UEqn.H"
                #include "YhEqn.H"

                // --- Pressure corrector loop
                while (pimple.correct())
                {
                    #include "pEqn.H"
                }

                if (pimple.turbCorr())
                {
                    turbulence->correct();
                }
            }

            rho = thermo.rho();
        }

        runTime.write();
        
        Info<< "\nExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
