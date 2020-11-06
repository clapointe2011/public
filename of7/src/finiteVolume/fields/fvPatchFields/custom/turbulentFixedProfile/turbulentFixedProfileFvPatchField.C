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

\*---------------------------------------------------------------------------*/

#include "turbulentFixedProfileFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::turbulentFixedProfileFvPatchField<Type>::turbulentFixedProfileFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    ranGen_(label(0)),
    fluctuationScale_(Zero),
    alpha_(0.1),
    curTimeIndex_(-1),
    profile_(),
    dir_(Zero),
    origin_(0)
{}


template<class Type>
Foam::turbulentFixedProfileFvPatchField<Type>::turbulentFixedProfileFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& fld
)
:
    fixedValueFvPatchField<Type>(p, iF, fld),
    ranGen_(label(0)),
    fluctuationScale_(Zero),
    alpha_(0.1),
    curTimeIndex_(-1),
    profile_(),
    dir_(Zero),
    origin_(0)
{}


template<class Type>
Foam::turbulentFixedProfileFvPatchField<Type>::turbulentFixedProfileFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict, false),
    ranGen_(label(0)),
    fluctuationScale_(pTraits<Type>(dict.lookup("fluctuationScale"))),
    alpha_(dict.lookupOrDefault<scalar>("alpha", 0.1)),
    curTimeIndex_(-1),
    profile_(Function1<Type>::New("profile", dict)),
    dir_(dict.lookup("direction")),
    origin_(readScalar(dict.lookup("origin")))
{
    if (mag(dir_) < small)
    {
        FatalErrorInFunction
            << "magnitude Direction must be greater than zero"
            << abort(FatalError);
    }

    // Ensure direction vector is normalized
    dir_ /= mag(dir_);

    // Evaluate profile
    this->evaluate();
}


template<class Type>
Foam::turbulentFixedProfileFvPatchField<Type>::turbulentFixedProfileFvPatchField
(
    const turbulentFixedProfileFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(p, iF),  // Don't map
    ranGen_(label(0)),
    fluctuationScale_(ptf.fluctuationScale_),
    alpha_(ptf.alpha_),
    curTimeIndex_(-1),
    profile_(ptf.profile_, false),
    dir_(ptf.dir_),
    origin_(ptf.origin_)
{
    // Evaluate profile since value not mapped
    this->evaluate();
}


template<class Type>
Foam::turbulentFixedProfileFvPatchField<Type>::turbulentFixedProfileFvPatchField
(
    const turbulentFixedProfileFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    ranGen_(ptf.ranGen_),
    fluctuationScale_(ptf.fluctuationScale_),
    alpha_(ptf.alpha_),
    curTimeIndex_(-1),
    profile_(ptf.profile_, false),
    dir_(ptf.dir_),
    origin_(ptf.origin_)
{}


template<class Type>
Foam::turbulentFixedProfileFvPatchField<Type>::turbulentFixedProfileFvPatchField
(
    const turbulentFixedProfileFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    ranGen_(ptf.ranGen_),
    fluctuationScale_(ptf.fluctuationScale_),
    alpha_(ptf.alpha_),
    curTimeIndex_(-1),
    profile_(ptf.profile_, false),
    dir_(ptf.dir_),
    origin_(ptf.origin_)
{
    // Evaluate the profile if defined
    if (ptf.profile_.valid())
    {
        this->evaluate();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::turbulentFixedProfileFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const scalarField dirCmpt((dir_ & this->patch().Cf()) - origin_);

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        Field<Type> patchField = *this;

        Field<Type> randomField(this->size());

        forAll(randomField, facei)
        {
            randomField[facei] = ranGen_.sample01<Type>();
        }

	fvPatchField<Type>::operator==
        (
            profile_->value(dirCmpt)
          + cmptMultiply
            (
                randomField - 0.5*pTraits<Type>::one,
                fluctuationScale_
            )*mag(profile_->value(dirCmpt))
        );

        // Correction-factor to compensate for the loss of RMS fluctuation
        // due to the temporal correlation introduced by the alpha parameter.
        /*scalar rmsCorr = sqrt(12*(2*alpha_ - sqr(alpha_)))/alpha_;

        patchField=
            (1 - alpha_)*patchField
          + alpha_*
            (
                profile_->value(dirCmpt)
              + rmsCorr*cmptMultiply
                (
                    randomField - 0.5*pTraits<Type>::one,
                    fluctuationScale_
                )*mag(profile_->value(dirCmpt))
	    );

        curTimeIndex_ = this->db().time().timeIndex();*/
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::turbulentFixedProfileFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    writeEntry(os, profile_());
    writeEntry(os, "fluctuationScale", fluctuationScale_);
    writeEntry(os, "direction", dir_);
    writeEntry(os, "origin", origin_);
    writeEntry(os, "value", *this);
}


// ************************************************************************* //
