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
    wilFireFoam

Author
    Adapted from original fireFoam code for use with AMR by Caelan Lapointe
    and Peter Hamlington.

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
#include "dynamicRefineBalanceMultiRegionFvMesh.H"
#include "cellSet.H"
#include "CorrectPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool run
(
    argList& argss,
    bool& balance,
    label& loopNumber,
    scalar& executionTime,
    scalar& clockTime
);

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    #include "setRootCase.H"

    label loopNumber = 1;
    scalar executionTime = 0.0;
    scalar clockTime = 0.0;
    bool balance = false;
    bool restart = false;
    do
    {
        restart = run
        (
            args,
            balance,
            loopNumber,
            executionTime,
            clockTime
        );
    } while(restart);

    Info<< "End" << endl;

    return 0;
}

bool run
(
    argList& args,
    bool& balance,
    label& loopNumber,
    scalar& executionTime,
    scalar& clockTime
)
{
    Info<< "\nStarting loop " << loopNumber << nl << endl;

    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createCustomControls.H"
    #include "createFieldRefs.H"
    
    turbulence->validate();

    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"
    
    #include "readPyrolysisTimeControls.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    
    scalar startRunIndex = runTime.timeIndex();
    scalar currentTime = runTime.time().value();
    scalar endTime = runTime.time().endTime().value();

    //Scaled Qdot field
    dimensionedScalar minsq = min(Qdot);
    dimensionedScalar maxsq = max(Qdot);
    volScalarField sQ("sQ", (Qdot - minsq)/(maxsq - minsq + dimensionedScalar(Qdot.dimensions(), small)));
    //Scaled grad sQ
    volScalarField vQdot = mag(fvc::grad(Qdot));
    dimensionedScalar minvsq = min(vQdot);
    dimensionedScalar maxvsq = max(vQdot);
    volScalarField varsQ("varsQ", (vQdot - minvsq)/(maxvsq - minvsq + dimensionedScalar(vQdot.dimensions(), small)));
    //Scaled enstrophy
    volScalarField enst = 0.5*magSqr(fvc::curl(U));
    dimensionedScalar minenst = min(enst);
    dimensionedScalar maxenst = max(enst);
    volScalarField sEnst("sEnst", (enst - minenst)/(maxenst - minenst + dimensionedScalar(enst.dimensions(), small)));

    while (runTime.run())
    {
        if (solvePrimaryRegion && (runTime.timeIndex() - startRunIndex) > startDyMIndex)
        {
            #include "updateMesh.H"
        }

        #include "calcCourantNo.H"
        #include "solidRegionDiffusionNo.H"
        #include "setMultiRegionDeltaT.H"
        #include "setDeltaT.H"

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
            if (!combInPimple) combustion->correct();
            if (!radInPimple) radiation->correct();

            #include "rhoEqn.H"

            // --- PIMPLE loop
            while (pimple.loop())
	    {
                #include "UEqn.H"
                
                if (!turbInPiso)
                {
                    #include "solveTurbulence.H"
                }

                if (!energyInPiso)
                {
                    #include "solveEnergy.H"
                }

                // --- Pressure corrector loop
                label pisoCorr = 1;
                while (pimple.correct())
                {
                    Info << "PISO: iteration " << pisoCorr << endl;

                    if (turbInPiso)
                    {
                        #include "solveTurbulence.H"
                    }

                    if (energyInPiso)
                    {
                        #include "solveEnergy.H"
                    }

                    #include "pEqn.H"
                    pisoCorr++;
                }
            }

	    rho = thermo.rho();
	}

        minsq = min(Qdot);
        maxsq = max(Qdot);
        sQ = (Qdot - minsq)/(maxsq - minsq + dimensionedScalar(Qdot.dimensions(), small));

        vQdot = mag(fvc::grad(Qdot));
        minvsq = min(vQdot);
        maxvsq = max(vQdot);
        varsQ = (vQdot - minvsq)/(maxvsq - minvsq + dimensionedScalar(vQdot.dimensions(), small));

        enst = 0.5*magSqr(fvc::curl(U));
        minenst = min(enst);
        maxenst = max(enst);
        sEnst = (enst - minenst)/(maxenst - minenst + dimensionedScalar(enst.dimensions(), small));

        runTime.write();
        
        Info<< "\nExecutionTime = " << executionTime + runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << clockTime + runTime.elapsedClockTime() << " s"
            << nl << endl;

        if (runTime.outputTime() && (currentTime != endTime))
	{
            balance = mesh.balance();
            if (balance)
            {
                executionTime += runTime.elapsedCpuTime();
                clockTime += runTime.elapsedClockTime();
                currentTime = runTime.time().value();
                runTime.setEndTime(currentTime);
            }
        }
    }

    if (balance)
    {
        mesh.redistribute();
        runTime.setEndTime(endTime);
        loopNumber++;
        balance = false;
        return true;
    }

    return false;
}



// ************************************************************************* //
