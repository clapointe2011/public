{
    volScalarField& rDeltaT = trDeltaT();

    volScalarField rDeltaT0("rDeltaT0", rDeltaT);

    const dictionary& pimpleDict = pimple.dict();

    // Maximum flow Courant number
    scalar maxCo(readScalar(pimpleDict.lookup("maxCo")));

    // Maximum time scale
    scalar maxDeltaT(pimpleDict.lookupOrDefault<scalar>("maxDeltaT", GREAT));

    // Smoothing parameter (0-1) when smoothing iterations > 0
    scalar rDeltaTSmoothingCoeff
    (
        pimpleDict.lookupOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.1)
    );

    // Damping coefficient (1-0)
    scalar rDeltaTDampingCoeff
    (
        pimpleDict.lookupOrDefault<scalar>("rDeltaTDampingCoeff", 1)
    );

    // Maximum change in cell temperature per iteration
    // (relative to previous value)
    scalar alphaTemp(pimpleDict.lookupOrDefault("alphaTemp", 0.05));

    if ( CourantType == "cellCourant" )
    {
        // Set the reciprocal time-step from the local Courant number

        rDeltaT.dimensionedInternalField() = max
        (
            1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
            (
                0.5 * fvc::surfaceSum(mag(phi))().dimensionedInternalField() 
                / mesh.V() / rho.dimensionedInternalField()
            ) / maxCo
        );
    }
    else if ( CourantType == "faceCourant" )
    {
        if (mesh.nInternalFaces())
        {
            
            surfaceScalarField phiOverRho =
            mag
            (
                phi / (p_pos * psi_pos * a_pos + p_neg * psi_neg * a_neg)
            );
            
            surfaceScalarField SfUfbyDelta =
                mesh.surfaceInterpolation::deltaCoeffs()*phiOverRho*2.0;
            
            forAll(mesh.V(), iCell)
            {
                const labelList& cellFaces = mesh.cells()[iCell];
                
                scalar maxCellRDeltaT = -1.0;
                scalar faceRDeltaT = maxCellRDeltaT;
                label faceId = -1;
                forAll(cellFaces, jFace)
                {
                    faceId = cellFaces[jFace];
                    if (mesh.isInternalFace(faceId))
                    {
                        faceRDeltaT = SfUfbyDelta.internalField()[faceId] / mesh.magSf().internalField()[faceId] / maxCo;
                    }
                    if (faceRDeltaT > maxCellRDeltaT)
                    {
                        maxCellRDeltaT = faceRDeltaT;
                    }
                }
                rDeltaT.internalField()[iCell] = maxCellRDeltaT;
            }
            forAll(SfUfbyDelta.boundaryField(), kPatch)
            {
                if(polyPatch::constraintType(mesh.boundaryMesh().types()[kPatch]))
                {
                }
                else
                {
                    const fvsPatchField<scalar>& pSfUfbyDelta = SfUfbyDelta.boundaryField()[kPatch];
                    const labelList& patchCells = mesh.boundary()[kPatch].faceCells();
                    scalar faceRDeltaT = -1;
                    forAll(pSfUfbyDelta, qFace)
                    {
                        faceRDeltaT = pSfUfbyDelta[qFace] / mesh.magSf().boundaryField()[kPatch][qFace] / maxCo;
                        if (faceRDeltaT > rDeltaT.internalField()[patchCells[qFace]])
                        {
                            rDeltaT.internalField()[patchCells[qFace]] = faceRDeltaT;
                        }
                    }
                }
            }
        }
        
    }
    else
    {
        FatalErrorIn
        (
            "setRDeltaT.H:"
        )
        << "Wrong type of Courant criterion: " << CourantType
        << endl << " must be one of:" 
        << endl << "1) cellCourant"
        << endl << "2) faceCourant"
        << endl << abort(FatalError);
    }

    // Update the boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    // Reaction source time scale
    if (alphaTemp < 1.0)
    {
        volScalarField::DimensionedInternalField rDeltaTT
        (
            mag(reaction->Sh())/(alphaTemp*rho*thermo.Cp()*T)
        );

        rDeltaT.dimensionedInternalField() = max
        (
            rDeltaT.dimensionedInternalField(),
            rDeltaTT
        );
    }

    if (rDeltaTSmoothingCoeff < 1.0)
    {
        fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
    }

    if
    (
        rDeltaTDampingCoeff < 1.0
     && runTime.timeIndex() > runTime.startTimeIndex() + 1
    )
    {
        rDeltaT = max
        (
            rDeltaT,
            (scalar(1.0) - rDeltaTDampingCoeff)*rDeltaT0
        );
    }

    rDeltaT.correctBoundaryConditions();

    Info<< "Flow time scale min/max = "
        << gMin(1/(rDeltaT.internalField() + 1.0e-12))
        << ", " << gMax(1/(rDeltaT.internalField() + 1.0e-12)) << endl;
}
