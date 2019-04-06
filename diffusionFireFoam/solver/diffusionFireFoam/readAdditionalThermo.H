/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
Global
    readAdditionalThermo
Description
\*---------------------------------------------------------------------------*/

const scalar lewisNo(thermo.lookupOrDefault<scalar>("lewisNo", 1));

if (lewisNo!=1) 
{
    Info<< "Lewis number: " << lewisNo  << "\n" << nl;
}

volScalarField dEff
(
    "dEff", turbulence->alphaEff()/(lewisNo + SMALL)
);

/*// moleculer mass diffusivity [m^2/s]
dimensionedScalar DM
(
    "DM",
     dimensionSet(0, 2, -1 , 0, 0),
     0
);
volScalarField dEff
(
    "dEff",
    turbulence->alphaEff()
);
const bool constDiffusivity(thermo.lookupOrDefault<bool>("constDiffusivity", false));
if (constDiffusivity)
{
    DM = dimensionedScalar(thermo.lookup("DM"));
    Info<< "Mass species diffusivity: " << DM  << "\n" << nl;
}
else
{
    if (lewisNo!=1) 
    {
        Info<< "Lewis number: " << lewisNo  << "\n" << nl;
    }
}*/

// ************************************************************************* //
