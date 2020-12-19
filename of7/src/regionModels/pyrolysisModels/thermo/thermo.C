/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2019 OpenCFD Ltd
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

\*---------------------------------------------------------------------------*/

#include "thermo.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "absorptionEmissionModel.H"
#include "fvm.H"
#include "fvcLaplacian.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace pyrolysisModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermo, 0);
addToRunTimeSelectionTable(pyrolysisModel, thermo, mesh);
addToRunTimeSelectionTable(pyrolysisModel, thermo, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void thermo::readThermoControls()
{
    const dictionary& solution = this->solution().subDict("SIMPLE");
    solution.lookup("nNonOrthCorr") >> nNonOrthCorr_;
    time().controlDict().lookup("maxDi") >> maxDiff_;
}


bool thermo::read()
{
    if (pyrolysisModel::read())
    {
        readThermoControls();
        return true;
    }

    return false;
}


bool thermo::read(const dictionary& dict)
{
    if (pyrolysisModel::read(dict))
    {
        readThermoControls();
        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermo::thermo
(
    const word& modelType,
    const fvMesh& mesh,
    const word& regionType
)
:
    pyrolysisModel(modelType, mesh, regionType),
    solidThermo_(solidThermo::New(regionMesh())),
    radiation_(radiationModel::New(solidThermo_->T())),
    h_(solidThermo_->he()),
    nNonOrthCorr_(-1),
    maxDiff_(10)
{
    if (active())
    {
         read();
    }
}


thermo::thermo
(
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& regionType
)
:
    pyrolysisModel(modelType, mesh, dict, regionType),
    solidThermo_(solidThermo::New(regionMesh())),
    radiation_(radiationModel::New(solidThermo_->T())),
    h_(solidThermo_->he()),
    nNonOrthCorr_(-1),
    maxDiff_(10)
{
    if (active_)
    {
        read(dict);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

thermo::~thermo()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar thermo::solidRegionDiffNo() const
{
    scalar DiNum = -great;

    if (regionMesh().nInternalFaces() > 0)
    {
        surfaceScalarField KrhoCpbyDelta
        (
            sqr(regionMesh().surfaceInterpolation::deltaCoeffs())
           *fvc::interpolate(kappa())
           /fvc::interpolate(Cp()*rho())
        );

        DiNum = max(KrhoCpbyDelta.primitiveField())*time().deltaTValue();
    }

    return DiNum;
}


scalar thermo::maxDiff() const
{
    return maxDiff_;
}

void thermo::preEvolveRegion()
{
    pyrolysisModel::preEvolveRegion();
}


void thermo::evolveRegion()
{
    Info<< "\nEvolving pyrolysis in region: " << regionMesh().name() << endl;

    for (int nonOrth=0; nonOrth<=nNonOrthCorr_; nonOrth++)
    {
        tmp<volScalarField> alpha(solidThermo_->alpha());
        fvScalarMatrix hEqn
        (
            fvm::ddt(rho(), h_)
          - fvm::laplacian(alpha, h_)
          + fvc::laplacian(alpha, h_)
          - fvc::laplacian(kappa(), T())
        );

        hEqn.relax();
        hEqn.solve();
    }

    solidThermo_->correct();

    Info<< "pyrolysis min/max(T) = "
        << gMin(solidThermo_->T().primitiveField())
        << ", "
        << gMax(solidThermo_->T().primitiveField())
        << nl << endl;
}


const volScalarField& thermo::rho() const
{
    return solidThermo_->rho();
}


const volScalarField& thermo::T() const
{
    return solidThermo_->T();
}


const tmp<volScalarField> thermo::Cp() const
{
    return solidThermo_->Cp();
}


tmp<volScalarField> thermo::kappaRad() const
{
    return radiation_->absorptionEmission().a();
}


tmp<volScalarField> thermo::kappa() const
{
    return solidThermo_->kappa();
}


const surfaceScalarField& thermo::phiGas() const
{
    FatalErrorInFunction
        << "phiGas field not available for " << type() << abort(FatalError);
    return surfaceScalarField::null();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //

