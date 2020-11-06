/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
\*---------------------------------------------------------------------------*/

#include "extendedEddyDissipation.H"
#include "uniformDimensionedFields.H"
#include "zeroGradientFvPatchFields.H"

namespace Foam
{
namespace combustionModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
extendedEddyDissipation<ReactionThermo, ThermoType>::extendedEddyDissipation
(
    const word& modelType,
    ReactionThermo& thermo,
    const compressibleTurbulenceModel& turb,
    const word& combustionProperties
)
:
    singleStepCombustion<ReactionThermo, ThermoType>
    (
        modelType,
        thermo,
        turb,
        combustionProperties
    ),
    C_(readScalar(this->coeffs().lookup("C_EDC"))),
    Cd_(readScalar(this->coeffs().lookup("C_Diff")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
extendedEddyDissipation<ReactionThermo, ThermoType>::~extendedEddyDissipation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
extendedEddyDissipation<ReactionThermo, ThermoType>::rtAdv() const
{
    return C_*this->turbulence().epsilon()/
        max(this->turbulence().k(),
        dimensionedScalar("SMALL",dimVelocity*dimVelocity,SMALL));
}

template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
extendedEddyDissipation<ReactionThermo, ThermoType>::rtDiff() const
{
    const volScalarField& YO2 = this->thermo().composition().Y("O2");
    const compressible::LESModel& lesModel =
        YO2.db().lookupObject<compressible::LESModel>
        (
            turbulenceModel::propertiesName
        );

    //Assume unity Prt and Lewis
    return Cd_*this->thermo().alpha()/this->rho()/sqr(lesModel.delta());
}

template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
extendedEddyDissipation<ReactionThermo, ThermoType>::rtBuoyAcc() const
{   
    const volScalarField& YO2 = this->thermo().composition().Y("O2");
    const compressible::LESModel& lesModel =
        YO2.db().lookupObject<compressible::LESModel>
        (
            turbulenceModel::propertiesName
        );
        
    const uniformDimensionedVectorField& g = YO2.db().lookupObject<uniformDimensionedVectorField>("g");

    return sqrt(mag(g)/(2*lesModel.delta()));
}

template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
extendedEddyDissipation<ReactionThermo, ThermoType>::rtChem() const
{
    tmp<volScalarField> ttChem
    (
        new volScalarField
        (
            IOobject
            (
                "ttChem",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("tmpChem",dimTime,Zero),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    volScalarField& tChem = ttChem.ref();
    tChem.ref() += dimensionedScalar("tmpChem",dimTime,1e-10);
    tChem.correctBoundaryConditions();
    
    return 1/tChem;
}

template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
extendedEddyDissipation<ReactionThermo, ThermoType>::rtFlame() const
{
    tmp<volScalarField> ttFlame
    (
        new volScalarField
        (
            IOobject
            (
                "ttFlame",
                this->mesh().time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("tmpFlame",dimTime,Zero),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    
    volScalarField& tFlame = ttFlame.ref();
    tFlame.ref() += dimensionedScalar("tmpFlame",dimTime,1e10);
    tFlame.correctBoundaryConditions();

    
    return 1/tFlame;
}

template<class ReactionThermo, class ThermoType>
void extendedEddyDissipation<ReactionThermo, ThermoType>::correct()
{
    this->wFuel_ ==
        dimensionedScalar("zero", dimMass/pow3(dimLength)/dimTime, 0.0);

    this->singleMixturePtr_->fresCorrect();

    const label fuelI = this->singleMixturePtr_->fuelIndex();

    const volScalarField& YFuel = this->thermo().composition().Y()[fuelI];

    const dimensionedScalar s = this->singleMixturePtr_->s();

    if (this->thermo().composition().contains("O2"))
    {
        const volScalarField& YO2 = this->thermo().composition().Y("O2");

        this->wFuel_ ==
              this->rho()
            * min(YFuel, YO2/s.value())
            * min( max( max( max(rtAdv(),rtDiff()), rtBuoyAcc()), rtFlame()), rtChem());
    }
}


template<class ReactionThermo, class ThermoType>
bool extendedEddyDissipation<ReactionThermo, ThermoType>::read()
{
    if (singleStepCombustion<ReactionThermo, ThermoType>::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

