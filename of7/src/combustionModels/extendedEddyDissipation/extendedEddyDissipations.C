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

#include "makeCombustionTypes.H"

#include "thermoPhysicsTypes.H"
#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"
#include "extendedEddyDissipation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Combustion models based on sensibleEnthalpy

makeCombustionTypesThermo
(
    extendedEddyDissipation,
    psiReactionThermo,
    gasHThermoPhysics
);

makeCombustionTypesThermo
(
    extendedEddyDissipation,
    psiReactionThermo,
    constGasHThermoPhysics
);

makeCombustionTypesThermo
(
    extendedEddyDissipation,
    rhoReactionThermo,
    gasHThermoPhysics
);

makeCombustionTypesThermo
(
    extendedEddyDissipation,
    rhoReactionThermo,
    constGasHThermoPhysics
);

// Combustion models based on sensibleInternalEnergy

makeCombustionTypesThermo
(
    extendedEddyDissipation,
    psiReactionThermo,
    gasEThermoPhysics
);

makeCombustionTypesThermo
(
    extendedEddyDissipation,
    psiReactionThermo,
    constGasEThermoPhysics
);

makeCombustionTypesThermo
(
    extendedEddyDissipation,
    rhoReactionThermo,
    gasEThermoPhysics
);

makeCombustionTypesThermo
(
    extendedEddyDissipation,
    rhoReactionThermo,
    constGasEThermoPhysics
);

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
