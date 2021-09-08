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
    rhoDiffusionFireFoam

Author
    Adapted from original fireFoam code for use with AMR by Caelan Lapointe
    and Peter Hamlington.

Description
    Transient solver for fires and turbulent diffusion flames.
    
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "radiationModel.H"
#include "rhoReactionThermo.H"
#include "CombustionModel.H"
#include "CorrectPhi.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "dynamicRefineFvMesh.H"
#include "CMULES.H"

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
    
    autoPtr<volVectorField> rhoU0;
    {
        rhoU0 = new volVectorField("rhoU0", rhoU);
    }
    
    autoPtr<volScalarField> divrhoU;
    if (correctPhi)
    {
        divrhoU = new volScalarField
        (
            "divrhoU",
            fvc::div(fvc::absolute(phi, rho, U))
        );
    }

    autoPtr<volScalarField> divrhoU0;
    if (correctPhi)
    {
        divrhoU0 = new volScalarField
        (
            "divrhoU0",
            divrhoU
        );
    }

    word ddtSchemeName(U.mesh().ddtScheme(U.name()));

    autoPtr<volScalarField> divrhoU00;
    autoPtr<volVectorField> rhoU00;
    if (correctPhi && (ddtSchemeName == "backward" || ddtSchemeName == "CrankNicolson"))
    {
        rhoU00 = new volVectorField("rhoU00", rhoU);

        divrhoU00 = new volScalarField
        (
            "divrhoU00",
            divrhoU
        );
    }

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
        
        if ((runTime.timeIndex() - startRunIndex) > startDyMIndex)
	{
            if (monitorConservation)
            {
                Info << "Mass : " << fvc::domainIntegrate(rho).value()<< endl;
                Info << "MomentumX : " << fvc::domainIntegrate(rho*U.component(0)).value() << endl;
                Info << "MomentumY : " << fvc::domainIntegrate(rho*U.component(1)).value() << endl;
                Info << "MomentumZ : " << fvc::domainIntegrate(rho*U.component(2)).value() << endl;
                Info << "Enthalpy : " << fvc::domainIntegrate(rho*thermo.he()).value() << endl;
                Info << "Count : " << mesh.C().size() << nl << endl;
	    }

            #include "updateMesh.H"
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

            if (useMULES)
            {
                #include "YEEqnMULES.H"
            }
            else
            {
                #include "YEEqn.H"
            }

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
