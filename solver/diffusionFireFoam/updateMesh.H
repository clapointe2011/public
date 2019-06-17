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
    
    // Make the fluxes absolute
    fvc::makeAbsolute(phi, rho, U);
    fvc::makeAbsolute(phi.oldTime(), rho.oldTime(), U.oldTime());
    if (ddtSchemeName == "backward")
    {
        fvc::makeAbsolute(phi.oldTime().oldTime(), rho.oldTime().oldTime(), U.oldTime().oldTime());
    }
    
    // Store momentum to set phi/phi_0 for introduced faces.
    rhoU() = rho*U;
    rhoU_0() = rho.oldTime()*U.oldTime();
    if (ddtSchemeName == "backward")
    {
        rhoU_0_0() = rho.oldTime().oldTime()*U.oldTime().oldTime();
    }
    
    mesh.update();

    if (mesh.changing())
    {                   
        //CL explicitly update gh/ghf
        gh = (g & mesh.C()) - ghRef;
        ghf = (g & mesh.Cf()) - ghRef;

        MRF.update();
                        
        if (correctPhi)
        {
            //CL explicitly update phi/phi_0/phi_0_0
            phi = mesh.Sf() & fvc::interpolate(rhoU());
            phi.oldTime() = mesh.Sf() & fvc::interpolate(rhoU_0());
            if (ddtSchemeName == "backward")
            {
                phi.oldTime().oldTime() = mesh.Sf() & fvc::interpolate(rhoU_0_0());
            }
        
            #include "correctPhi.H"
                            
            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phi, rho, U);
            fvc::makeRelative(phi.oldTime(), rho.oldTime(), U.oldTime());
            if (ddtSchemeName == "backward")
            {
                fvc::makeRelative(phi.oldTime().oldTime(), rho.oldTime().oldTime(), U.oldTime().oldTime());
            }
        }
    }

    Info<< "\nExecution time for mesh.update() = "
        << runTime.elapsedCpuTime() - timeBeforeMeshUpdate
        << " s" << nl << endl;
}
