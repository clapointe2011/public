/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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
    Foam::turbulentFixedProfileFvPatchField

Description
    This boundary condition provides a fixed value profile condition.

Usage
    \table
        Property     | Description              | Required | Default value
        profile      | Profile Function1        | yes |
        direction    | Profile direction        | yes |
        origin       | Profile origin           | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            turbulentFixedProfile;
        profile         csvFile;

        profileCoeffs
        {
            nHeaderLine         0;          // Number of header lines
            refColumn           0;          // Reference column index
            componentColumns    (1 2 3);    // Component column indices
            separator           ",";        // Optional (defaults to ",")
            mergeSeparators     no;         // Merge multiple separators
            file                "Uprofile.csv";  // name of csv data file
            outOfBounds         clamp;      // Optional out-of-bounds handling
            interpolationScheme linear;     // Optional interpolation scheme
        }
        direction        (0 1 0);
        origin           0;
    }
    \endverbatim

    Example setting a parabolic inlet profile for the PitzDaily case:
    \verbatim
    inlet
    {
        type            turbulentFixedProfile;

        profile         polynomial
        (
            ((1 0 0)        (0 0 0))
            ((-6200 0 0)    (2 0 0))
        );
        direction       (0 1 0);
        origin          0.0127;
    }
    \endverbatim

Note
    The profile entry is a Function1 type.  The example above gives the
    usage for supplying csv file.

See also
    Foam::fixedValueFvPatchField
    Foam::Function1Types
    Foam::timeVaryingMappedFixedValueFvPatchField

SourceFiles
    turbulentFixedProfileFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentFixedProfileFvPatchField_H
#define turbulentFixedProfileFvPatchField_H

#include "Random.H"
#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class turbulentFixedProfileFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class turbulentFixedProfileFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private Data

        //- Random number generator
        Random ranGen_;

        //- Fluctuation scale
        Type fluctuationScale_;

        //- Fraction of RMS component to apply to last time step values
        scalar alpha_;

        //- Current time index (used for updating)
        label curTimeIndex_;

        //- Profile data
        autoPtr<Function1<Type>> profile_;

        //- Profile direction
        vector dir_;

        //- Profile origin
        scalar origin_;


public:

    //- Runtime type information
    TypeName("turbulentFixedProfile");


    // Constructors

        //- Construct from patch and internal field
        turbulentFixedProfileFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch and internal field and patch field
        turbulentFixedProfileFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const Field<Type>& fld
        );

        //- Construct from patch, internal field and dictionary
        turbulentFixedProfileFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given turbulentFixedProfileFvPatchField
        //  onto a new patch
        turbulentFixedProfileFvPatchField
        (
            const turbulentFixedProfileFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        turbulentFixedProfileFvPatchField
        (
            const turbulentFixedProfileFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new turbulentFixedProfileFvPatchField<Type>(*this)
            );
        }

        //- Copy constructor setting internal field reference
        turbulentFixedProfileFvPatchField
        (
            const turbulentFixedProfileFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new turbulentFixedProfileFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "turbulentFixedProfileFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
