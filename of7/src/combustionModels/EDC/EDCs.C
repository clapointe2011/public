/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"
#include "EDC.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

template<>
const char* Foam::NamedEnum
<
    Foam::combustionModels::EDCversions,
    4
>::names[] =
{
    "v1981",
    "v1996",
    "v2005",
    "v2016"
};

const Foam::NamedEnum<Foam::combustionModels::EDCversions, 4>
    Foam::combustionModels::EDCversionNames;

const Foam::combustionModels::EDCversions
Foam::combustionModels::EDCdefaultVersion
(
    Foam::combustionModels::EDCversions::v2005
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makeCombustionTypes(EDC, psiReactionThermo);
makeCombustionTypes(EDC, rhoReactionThermo);

}

// ************************************************************************* //
