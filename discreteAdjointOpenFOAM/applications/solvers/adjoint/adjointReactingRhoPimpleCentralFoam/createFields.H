Info<< "Creating reaction model\n" << endl;

autoPtr<combustionModels::psiCombustionModel> reaction
(
    combustionModels::psiCombustionModel::New(mesh)
);

psiReactionThermo& thermo = reaction->thermo();
thermo.validate(args.executable(), "h", "e");
thermo.correct();

volScalarField& h = thermo.he();
volScalarField& p = thermo.p();

volScalarField& psi = const_cast<volScalarField&>(thermo.psi());
const volScalarField& T = thermo.T();

basicMultiComponentMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

word inertSpecie(thermo.lookup("inertSpecie"));

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField K
(
    "K",
    0.5*magSqr(U)
);

//CL create fields for gamma, Cp, Cv
volScalarField Cp
(
    "Cp",
    thermo.Cp()
);

volScalarField Cv
(
    "Cv",
    thermo.Cv()
);

volScalarField gamma
(
    "gamma",
    thermo.Cp()/thermo.Cv()
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

volScalarField rhoHat
(
    IOobject
    (
        "rhoHat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho - psi*p
);

dimensionedScalar rhoMax
(
    dimensionedScalar::lookupOrDefault
    (
        "rhoMax",
        pimple.dict(),
        dimDensity,
        GREAT
    )
);

dimensionedScalar rhoMin
(
    dimensionedScalar::lookupOrDefault
    (
        "rhoMin",
        pimple.dict(),
        dimDensity,
        SMALL
    )
);

dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);

volScalarField dpdt
(
    "dpdt",
    p / runTime.deltaT() * 0.0
);

volScalarField dissip
(
    "dissip",
    dpdt
);

volVectorField gradp
(
    "gradp",
    fvc::grad(p)
);

volVectorField limGradP 
(
    "limGradP",
    fvc::grad(p)
);

DynamicList<label> badQualityCells(0);

surfaceScalarField pos
(
    IOobject
    (
        "pos",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("pos", dimless, 1.0)
);

surfaceScalarField neg
(
    IOobject
    (
        "neg",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("neg", dimless, -1.0)
);

surfaceScalarField phi
(
    IOobject
    (
	"phi",
	runTime.timeName(),
	mesh,
	IOobject::READ_IF_PRESENT,
	IOobject::AUTO_WRITE
    ),
    linearInterpolate(rho*U) & mesh.Sf()
);

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

// Set the turbulence into the combustion model
reaction->setTurbulence(turbulence());

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll (Y,i)
{
    fields.add(Y[i]);
}

fields.add(thermo.he());

volScalarField dQ
(
    IOobject
    (
        "dQ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("dQ", dimEnergy/dimTime, 0.0)
);

#include "createOptFields.H"

