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
    
    Adapted from fireFoam-dev code for use with AMR by Caelan Lapointe

\*---------------------------------------------------------------------------*/

#include "eddyDissipationModel.H"

namespace Foam
{
namespace combustionModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
eddyDissipationModel<ReactionThermo, ThermoType>::eddyDissipationModel
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
    Cd_(readScalar(this->coeffs().lookup("C_Diff"))),
    Cstiff_(readScalar(this->coeffs().lookup("C_Stiff")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
eddyDissipationModel<ReactionThermo, ThermoType>::~eddyDissipationModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
eddyDissipationModel<ReactionThermo, ThermoType>::rtTurb() const
{
    return C_*this->turbulence().epsilon()/
              max(this->turbulence().k(),
              dimensionedScalar("SMALL",dimVelocity*dimVelocity,SMALL));
}

template<class ReactionThermo, class ThermoType>
Foam::tmp<Foam::volScalarField>
eddyDissipationModel<ReactionThermo, ThermoType>::rtDiff() const
{
    const volScalarField& YO2 = this->thermo().composition().Y("O2");
    const compressible::LESModel& lesModel =
        YO2.db().lookupObject<compressible::LESModel>
        (
            turbulenceModel::propertiesName
        );

    return Cd_*this->thermo().alpha()/this->rho()/sqr(lesModel.delta());
}

template<class ReactionThermo, class ThermoType>
void eddyDissipationModel<ReactionThermo, ThermoType>::correct()
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
        
        volScalarField rt(max(rtTurb(),rtDiff())); 

        this->wFuel_ ==
              this->rho()
            * min(YFuel, YO2/s.value())
            / this->mesh().time().deltaT() / Cstiff_
            * (1 - exp(- Cstiff_*this->mesh().time().deltaT() * rt));
    }
}


template<class ReactionThermo, class ThermoType>
bool eddyDissipationModel<ReactionThermo, ThermoType>::read()
{
    if (singleStepCombustion<ReactionThermo, ThermoType>::read())
    {
        this->coeffs().lookup("C_EDC") >> C_ ;
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
