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

Author
    Adapted from original fireFoam code for use with AMR by Caelan Lapointe
    and Peter Hamlington.

Description
    Transient solver for fires and turbulent diffusion flames with reacting
    particle clouds, surface film and pyrolysis modelling.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "turbulentFluidThermoModel.H"
#include "radiationModel.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "CorrectPhi.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "motionSolver.H"
#include "dynamicRefineFvMesh.H"
#include "bound.H"
#include "backwardDdtScheme.H"

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
    
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createRDeltaT.H"
    #include "createCustomControls.H"
    
    autoPtr<volVectorField> rhoU;
    {
        rhoU = new volVectorField("rhoU", rho*U);
    }
    
    autoPtr<volVectorField> rhoU_0;
    {
        rhoU_0 = new volVectorField("rhoU_0", rhoU);
    }
    
    autoPtr<volVectorField> rhoU_0_0;
    {
        rhoU_0_0 = new volVectorField("rhoU_0_0", rhoU);
    }
    
    word ddtSchemeName(U.mesh().ddtScheme(U.name()));

    turbulence->validate();

    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

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
        if (correctPhi && ddtSchemeName == "backward")
        {
            divrhoU_0_0 = new volScalarField
            (
                "divrhoU_0_0",
                fvc::div(fvc::absolute(phi.oldTime().oldTime(), rho.oldTime().oldTime(), U.oldTime().oldTime()))
            );
        }
        
        if ((runTime.timeIndex() - startRunIndex) > startDyMIndex)
        {
            #include "updateMesh.H"
        }
        
        autoPtr<surfaceScalarField> UBlendingFactor;
        if (localBlendedU)
        {
            volScalarField level
            (
                IOobject
                (
                    "level",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("tmp", dimless, 1.0),
                zeroGradientFvPatchScalarField::typeName
            );

            labelList cellLevel = mesh.lookupObject<labelIOList>("cellLevel");

            forAll(level, celli)
            {
                level[celli] = cellLevel[celli];
            }

            level /= max(level);

            UBlendingFactor = new surfaceScalarField
            (
                "UBlendingFactor",
                fvc::interpolate(level)
            );
        }
        
        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;
        
        Info<< "Time = " << runTime.timeName() << nl << endl;
        
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

        runTime.write();

        Info<< "\nExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
