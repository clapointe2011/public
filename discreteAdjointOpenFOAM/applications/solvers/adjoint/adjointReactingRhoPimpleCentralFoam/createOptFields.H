dictionary optDict = mesh.solutionDict().subDict("OPTIMIZATION");

autoPtr<wordList> costFunctionPatches;
costFunctionPatches.reset(new wordList(0));

if (optDict.found("costFunctionPatches"))
{
    costFunctionPatches.reset
    (
       new wordList(optDict.lookup("costFunctionPatches"))
    );
}else{
    Info << "Warning: Keyword costFunctionPatches not found in fvSolution/OPTIMIZATION" << endl;
    Info << "Example:\nOPTIMIZATION{\n  costFunctionPatches (inlet outlet);\n}" << endl;
}

// Towara: create field to hold resistance term alpha
dimensionedScalar zeroSensitivity("0", dimVelocity*dimVelocity, 0.0);
dimensionedScalar zeroAlpha("0", dimless/dimTime, 0.0);

//const labelList& inletCells = mesh.boundary()["inlet"].faceCells(); //CL
//const labelList& outletCells = mesh.boundary()["outlet"].faceCells(); //CL

//CL scalars for optimization
scalar relTol = readScalar(optDict.lookup("relTol")); //CL
scalar tapedIter = readScalar(optDict.lookup("tapedIter")); //CL
scalar optEndIter = readScalar(optDict.lookup("optEndIter")); //CL
scalar initCutoff = readScalar(optDict.lookup("initCutoff")); //CL
//scalar initCutoff1 = readScalar(optDict.lookup("initCutoff1")); //CL
//scalar initCutoff2 = readScalar(optDict.lookup("initCutoff2")); //CL
scalar finalCutoff = readScalar(optDict.lookup("finalCutoff")); //CL

dimensionedScalar lambda(optDict.lookup("lambda")); //CL
dimensionedScalar alphaMax(optDict.lookup("alphaMax")); //CL

bool passiveInitialization = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<bool>("passiveInitialization",false); //CL
bool passiveFinalization = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<bool>("passiveFinalization",false); //CL
bool cumulativeSolution = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<bool>("cumulativeSolution",false); //CL
//bool fromSteady = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<bool>("fromSteady",false); //CL
bool writeOptFields = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<bool>("writeOptFields",false); //CL
bool buildFrom = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<bool>("buildFrom",false); //CL
bool scaleByVol = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<bool>("scaleByVol",false); //CL
bool limitMinAlpha = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<bool>("limitMinAlpha",false); //CL
bool binaryResField = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<bool>("binaryAlpha",false); //CL
//bool proximityDist = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<bool>("proximityDist",false); //CL
bool proximityNeighbor = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<bool>("proximityNeighbor",false); //CL
bool smoothResField = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<bool>("smoothAlpha",false); //CL
bool fillGaps = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<bool>("fillGaps",false); //CL
bool frozenTurbulence = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<bool>("frozenTurbulence",false);

//scalar maxSeparation = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("maxSeparation",1); //CL
int writeOptFreq = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<int>("writeOptFreq",1); //CL
int smoothFreq = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<int>("smoothFreq",1); //CL
scalar MAfactor = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("MAfactor",0); //CL
scalar SAfactor = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("SAfactor",0.5); //CL
//scalar PSfactor = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("PSfactor",0); //CL
scalar SAweight = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("SAweight",0.5); //CL
scalar PRweight = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("PRweight",1); //CL

scalar tol = 0;
scalar delta = 0;
scalar prevJ = 0;
scalar prevAlpha = 0;

//CL define and initialize sensitivity field
volScalarField sens
(
    IOobject
    (
        "sens",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    zeroSensitivity,
    zeroGradientFvPatchScalarField::typeName
);

//CL define and initialize resistance field
volScalarField alpha
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    //lambda*max(-sens, zeroSensitivity)

    mesh,
    zeroAlpha,
    zeroGradientFvPatchScalarField::typeName
);

volScalarField smoothAlpha
(
    "smoothAlpha",
    alpha
);
