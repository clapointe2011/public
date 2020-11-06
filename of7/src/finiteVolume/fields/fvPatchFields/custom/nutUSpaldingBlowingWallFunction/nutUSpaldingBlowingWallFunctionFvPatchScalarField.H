/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
Class
    Foam::nutUSpaldingBlowingWallFunctionFvPatchScalarField
Group
    grpWallFunctions
Description
    This boundary condition provides a turbulent kinematic viscosity condition
    when using wall functions for rough walls, based on velocity,  using
    Spalding's law to give a continuous nut profile to the wall (y+ = 0)
        \f[
            y^+ = u^+ + \frac{1}{E} \left[exp(\kappa u^+) - 1 - \kappa u^+\,
                - 0.5 (\kappa u^+)^2 - \frac{1}{6} (\kappa u^+)^3\right]
        \f]
    where
    \vartable
        y^+     | non-dimensional position
        u^+     | non-dimensional velocity
        \kappa  | Von Karman constant
    \endvartable
    \heading Patch usage
    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            nutUSpaldingBlowingWallFunction;
    }
    \endverbatim
SeeAlso
    Foam::nutWallFunctionFvPatchScalarField
SourceFiles
    nutUSpaldingBlowingWallFunctionFvPatchScalarField.C
\*---------------------------------------------------------------------------*/

#ifndef nutUSpaldingBlowingWallFunctionFvPatchScalarField_H
#define nutUSpaldingBlowingWallFunctionFvPatchScalarField_H

#include "nutUSpaldingWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class nutUSpaldingBlowingWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutUSpaldingBlowingWallFunctionFvPatchScalarField
:
    public nutUSpaldingWallFunctionFvPatchScalarField
{

     // Private data

        //- h over Cp (default = 10), unit (W/m^2/K)/(J/g/K)=g/s/m^2
        scalar hOCp_;

        //- ratio between phi_pyrolysate and phi_fuel (phiGas/phi), default = 2.5 for corrugated over C3H8
        scalar rPhi_;

protected:

    // Protected Member Functions

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNut() const;

public:

    //- Runtime type information
    TypeName("nutUSpaldingBlowingWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutUSpaldingBlowingWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutUSpaldingBlowingWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nutUSpaldingBlowingWallFunctionFvPatchScalarField
        //  onto a new patch
        nutUSpaldingBlowingWallFunctionFvPatchScalarField
        (
            const nutUSpaldingBlowingWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutUSpaldingBlowingWallFunctionFvPatchScalarField
        (
            const nutUSpaldingBlowingWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nutUSpaldingBlowingWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nutUSpaldingBlowingWallFunctionFvPatchScalarField
        (
            const nutUSpaldingBlowingWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new nutUSpaldingBlowingWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
