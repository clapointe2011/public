const volScalarField& psi = thermo.psi();
const volScalarField& T = thermo.T();
regionModels::surfaceFilmModel& surfaceFilm = tsurfaceFilm();
const label inertIndex(composition.species()[inertSpecie]);

autoPtr<volVectorField> rhoU;
rhoU = new volVectorField
(
    "rhoU", rho*U
);

autoPtr<volVectorField> rhoU0;
rhoU0 = new volVectorField
(
    "rhoU0", rhoU()
);

autoPtr<volVectorField> rhoU00;
rhoU00 = new volVectorField
(
    "rhoU00", rhoU()
);

autoPtr<volScalarField> divrhoU;
if (correctPhi)
{
    divrhoU = new volScalarField
    (
        "divrhoU",
        fvc::div(phi)
    );
}

autoPtr<volScalarField> divrhoU0;
if (correctPhi)
{
    divrhoU0 = new volScalarField
    (
        "divrhoU0",
        divrhoU()
    );
}

autoPtr<volScalarField> divrhoU00;
if (correctPhi)
{
    divrhoU00 = new volScalarField
    (
        "divrhoU00",
        divrhoU()
    );
}
