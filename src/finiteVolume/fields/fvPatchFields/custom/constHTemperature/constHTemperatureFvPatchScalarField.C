/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "constHTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOobjectList.H"
#include "solidThermo.H"

#define DEBUG(x) std::cout << "["<< __FILE__ << ":" << __LINE__ << "] "<< #x " = " << x << std::endl;
#define TRACE(s) std::cout << "["<< __FILE__ << ":" << __LINE__ << "] "<< #s << std::endl; s;
#define DEBUGP(x) std::cout << "[p"<<Pstream::myProcNo()<<":"<< __FILE__ << ":" << __LINE__ << "] "<< #x " = " << x << std::endl;
#define TRACEP(s) std::cout << "[p"<<Pstream::myProcNo()<<":"<< __FILE__ << ":" << __LINE__ << "] "<< #s << std::endl; s;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::constHTemperatureFvPatchScalarField::
constHTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    Tinf_(p.size(), 0.0),
    h_(p.size(), 0.0)
{
    refValue() = 295;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::constHTemperatureFvPatchScalarField::
constHTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF),
    Tinf_("Tinf", dict, p.size()),
    h_("h", dict, p.size())
{
    refValue() = Tinf_;
    refGrad() = 0.0;
    valueFraction() = 0.0;

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            Field<scalar>("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(refValue());
    }
}

Foam::constHTemperatureFvPatchScalarField::
constHTemperatureFvPatchScalarField
(
    const constHTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper),
    Tinf_(ptf.Tinf_, mapper),
    h_(ptf.h_, mapper)
{}


Foam::constHTemperatureFvPatchScalarField::
constHTemperatureFvPatchScalarField
(
    const constHTemperatureFvPatchScalarField& tppsf
)
:
    mixedFvPatchField<scalar>(tppsf),
    Tinf_(tppsf.Tinf_),
    h_(tppsf.h_)
{}

Foam::constHTemperatureFvPatchScalarField::
constHTemperatureFvPatchScalarField
(
    const constHTemperatureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(tppsf, iF),
    Tinf_(tppsf.Tinf_),
    h_(tppsf.h_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constHTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    scalarField::autoMap(m);
    Tinf_.autoMap(m);
    h_.autoMap(m);
}


void Foam::constHTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchField<scalar>::rmap(ptf, addr);

    const constHTemperatureFvPatchScalarField& tiptf =
         refCast<const constHTemperatureFvPatchScalarField>(ptf);

    Tinf_.rmap(tiptf.Tinf_, addr);
    h_.rmap(tiptf.h_, addr);
}

void Foam::constHTemperatureFvPatchScalarField::updateCoeffs()
{

    if (this->updated())
    {
        return;
    }

//    scalarField K_ = patch().lookupPatchField<volScalarField, scalar>("K");
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const solidThermo& thermo =
        mesh.lookupObject<solidThermo>("thermophysicalProperties");

    scalarField K_(thermo.kappa(patch().index()));




    refValue() = Tinf_;
    refGrad() = 0.0;
    valueFraction() =
        1.0/(1.0 + K_/max(h_,SMALL)*patch().deltaCoeffs());

    mixedFvPatchField<scalar>::updateCoeffs();
}

// #################
void Foam::constHTemperatureFvPatchScalarField::
setTInf(const scalarField& TInf)
{


    /*DEBUGP(TInf.size());*/
    /*DEBUGP(Tinf_.size());*/
    if (TInf.size() != Tinf_.size())
    {
        FatalErrorIn
        (
            "Foam::constHTemperatureFvPatchScalarField"
            "setTInf"
            "("
                "const scalarField& "
            ")"
        )
            << " Patch size does not match in the setTinf function \n"
            << nl << nl
            << "There seems to be an issue in the BC specifications or in the code. "
            << " Please contact the code development team." << exit(FatalError);
    }
    Tinf_ = TInf;
}
// #####################

void Foam::constHTemperatureFvPatchScalarField::
write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    Tinf_.writeEntry("Tinf", os);
    h_.writeEntry("h", os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        constHTemperatureFvPatchScalarField
    );

}

// ************************************************************************* //
