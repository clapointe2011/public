
{
    scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();
    
    parcels.storeGlobalPositions();
                        
    // Test : disable refinement for some cells
    PackedBoolList& protectedCell = refCast<myDynamicRefineFvMesh>(mesh).protectedCell();

    if (protectedCell.empty())
    {
        protectedCell.setSize(mesh.nCells());
        protectedCell = 0;
    }
                        
    const fvPatchList& patches = mesh.boundary();
                        
    forAll(patches, patchi)
    {
	    const fvPatch& currPatch = patches[patchi];
	                        
	    if (protectMapped && mesh.boundary()[patchi].type() == "mappedWall")
	    {
	                        
	        forAll(currPatch,facei)
	        {
		        label faceCelli = currPatch.faceCells()[facei];
		        protectedCell[faceCelli] = 1;
	        }
	    }
    }
    
    if (addOnlyNearMapped)
    {
        cellSet protectedCells
        (
            IOobject
            (
                "protectedCells",
                "constant/polyMesh/sets",
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        
        forAll (mesh.C(), celli)
        {
            if (protectedCells[celli])
            {
                protectedCell[celli] = 1;
            }
        }
    }

    // Store flux divergences
    if (correctPhi)
    {
        divrhoU() = fvc::div(fvc::absolute(phi, rho, U));
        divrhoU_0() = fvc::div(fvc::absolute(phi.oldTime(), rho.oldTime(), U.oldTime()));
        if (ddtSchemeName != "Euler")
        {
            divrhoU_0_0() = fvc::div
            (
                fvc::absolute
                (
                    phi.oldTime().oldTime(),
                    rho.oldTime().oldTime(),
                    U.oldTime().oldTime()
                )
            );
        }
    }

    // Make the fluxes absolute
    fvc::makeAbsolute(phi, rho, U);
    fvc::makeAbsolute(phi.oldTime(), rho.oldTime(), U.oldTime());
    if (ddtSchemeName != "Euler")
    {
        fvc::makeAbsolute(phi.oldTime().oldTime(), rho.oldTime().oldTime(), U.oldTime().oldTime());
    }
    
    // Store momentum to set phi/phi_0 for introduced faces.
    rhoU() = rho*U;
    rhoU_0() = rho.oldTime()*U.oldTime();
    if (ddtSchemeName != "Euler")
    {
        rhoU_0_0() = rho.oldTime().oldTime()*U.oldTime().oldTime();
    }
    
    mesh.update();

    if (mesh.changing())
    {                   
        //CL explicitly update gh/ghf
        gh = (g & mesh.C()) - ghRef;
        ghf = (g & mesh.Cf()) - ghRef;
        
        //CL update ph_rgh
        /*if (pimple.dict().lookupOrDefault<bool>("hydrostaticInitialization", false))
        {
            surfaceScalarField ph_rghf = fvc::interpolate(ph_rgh);
            ph_rgh.clear();
            ph_rgh = fvc::average(ph_rghf);
            ph_rghf.clear();
        }*/

        MRF.update();
                        
        if (correctPhi)
        {
            //CL explicitly update phi/phi_0/phi_0_0
            phi = mesh.Sf() & fvc::interpolate(rhoU());
            phi.oldTime() = mesh.Sf() & fvc::interpolate(rhoU_0());
            if (ddtSchemeName != "Euler")
            {
                phi.oldTime().oldTime() = mesh.Sf() & fvc::interpolate(rhoU_0_0());
            }
        
            #include "correctPhi.H"
                            
            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phi, rho, U);
            fvc::makeRelative(phi.oldTime(), rho.oldTime(), U.oldTime());
            if (ddtSchemeName != "Euler")
            {
                fvc::makeRelative(phi.oldTime().oldTime(), rho.oldTime().oldTime(), U.oldTime().oldTime());
            }
        }
    }

    Info<< "\nExecution time for mesh.update() = "
        << runTime.elapsedCpuTime() - timeBeforeMeshUpdate
        << " s" << nl << endl;
}
