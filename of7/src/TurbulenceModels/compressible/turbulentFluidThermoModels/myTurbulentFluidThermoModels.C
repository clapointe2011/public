/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
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

#include "turbulentFluidThermoModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// -------------------------------------------------------------------------- //
// LES models
// -------------------------------------------------------------------------- //

#include "myDynamicKEqn.H"
makeLESModel(myDynamicKEqn);

#include "FLAMES_kEqn.H"
makeLESModel(FLAMES_kEqn);

#include "FLAMES_WHALE.H"
makeLESModel(FLAMES_WHALE);

#include "FLAMES_buoyantWHALE.H"
makeLESModel(FLAMES_buoyantWHALE);

#include "FLAMES_buoyantKEqn.H"
makeLESModel(FLAMES_buoyantKEqn);

#include "WHALE.H"
makeLESModel(WHALE);

#include "InagiKEqn.H"
makeLESModel(InagiKEqn);

#include "buoyantKEqn.H"
makeLESModel(buoyantKEqn);

#include "buoyantWHALE.H"
makeLESModel(buoyantWHALE);

// ************************************************************************* //
