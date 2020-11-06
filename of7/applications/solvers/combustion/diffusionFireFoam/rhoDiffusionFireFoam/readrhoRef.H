Info<< "\nReading rhoRef" << endl;

uniformDimensionedScalarField rhoRef
(
    IOobject
    (
        "rhoRef",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    dimensionedScalar(dimDensity, 0)
);
