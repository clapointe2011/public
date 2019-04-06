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
    #include "createDyMControls.H"
    
    autoPtr<volScalarField> rAU;
    {
        rAU = new volScalarField
        (
            IOobject
            (
                "rAU",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("rAU", dimTime/dimDensity, 1)
        );
    }
    
    autoPtr<volVectorField> rhoU;
    {
        rhoU = new volVectorField("rhoU", rho*U);
    }
    
    autoPtr<volVectorField> rhoU_0;
    {
        rhoU_0 = new volVectorField("rhoU_0", rhoU);
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
        #include "readDyMControls.H"
    
        autoPtr<volScalarField> divrhoU;
        if (correctPhi)
        {
            divrhoU = new volScalarField
            (
                "divrhoU",
                fvc::div(fvc::absolute(phi, rho, U))
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
            if (pimple.firstIter())
            {
                if ((runTime.timeIndex() - startRunIndex) > startDyMIndex)
                {
                    scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();
                    
                    if (protectOutlet)
                    {
                        // Test : disable refinement for some cells
                        PackedBoolList& protectedCell = refCast<dynamicRefineFvMesh>(mesh).protectedCell();

                        if (protectedCell.empty())
                        {
                            protectedCell.setSize(mesh.nCells());
                            protectedCell = 0;
                        }
                            
                        label patchLabel = mesh.boundaryMesh().findPatchID("outlet");
                        const fvPatch& patch = mesh.boundary()[patchLabel];
                           
                        forAll(patch,facei)
	                    {
	                        label faceCelli = patch.faceCells()[facei];
		                    protectedCell[faceCelli] = 1;
	                    }
	                }
                    
                    // Store momentum to set rhoUf for introduced faces.
                    rhoU() = rho*U;
                    rhoU_0() = rho.oldTime()*U.oldTime();

                    //CL use rhoU/rhoU0 to update phi/phi_0 as part of mesh.update()?
                    mesh.update();

                    if (mesh.changing())
                    {
                        phi = mesh.Sf() & fvc::interpolate(rhoU());
                        phi.oldTime() = mesh.Sf() & fvc::interpolate(rhoU_0());
                        
                        gh = (g & mesh.C()) - ghRef;
                        ghf = (g & mesh.Cf()) - ghRef;

                        MRF.update();
                        
                        if (correctPhi)
                        {
                            #include "correctPhi.H"
                            
                            // Make the fluxes relative to the mesh motion
                            fvc::makeRelative(phi, rho, U);
                        }
                    }

                    Info<< "\nExecution time for mesh.update() = "
                        << runTime.elapsedCpuTime() - timeBeforeMeshUpdate
                        << " s" << nl << endl;
                }
            }
        
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
