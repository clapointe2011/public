/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "fixedEnthalpyFluxTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
//#include "turbulenceModel.H"
#include "LESModel.H"
#include "IOobjectList.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedEnthalpyFluxTemperatureFvPatchScalarField::
fixedEnthalpyFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    phiName_("phi"),
    rhoName_("none"),
    Tinf_(p.size(), 0.0)
{
    refValue() = 295;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::fixedEnthalpyFluxTemperatureFvPatchScalarField::
fixedEnthalpyFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "none")),
    Tinf_("Tinf", dict, p.size())
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

Foam::fixedEnthalpyFluxTemperatureFvPatchScalarField::
fixedEnthalpyFluxTemperatureFvPatchScalarField
(
    const fixedEnthalpyFluxTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    Tinf_(mapper(ptf.Tinf_))
{}


Foam::fixedEnthalpyFluxTemperatureFvPatchScalarField::
fixedEnthalpyFluxTemperatureFvPatchScalarField
(
    const fixedEnthalpyFluxTemperatureFvPatchScalarField& tppsf
)
:
    mixedFvPatchField<scalar>(tppsf),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    Tinf_(tppsf.Tinf_)
{}

Foam::fixedEnthalpyFluxTemperatureFvPatchScalarField::
fixedEnthalpyFluxTemperatureFvPatchScalarField
(
    const fixedEnthalpyFluxTemperatureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(tppsf, iF),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    Tinf_(tppsf.Tinf_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedEnthalpyFluxTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchField::autoMap(m);
    m(Tinf_);
}


void Foam::fixedEnthalpyFluxTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchField<scalar>::rmap(ptf, addr);

    const fixedEnthalpyFluxTemperatureFvPatchScalarField& tiptf =
         refCast<const fixedEnthalpyFluxTemperatureFvPatchScalarField>(ptf);

    Tinf_.rmap(tiptf.Tinf_, addr);
}

void Foam::fixedEnthalpyFluxTemperatureFvPatchScalarField::updateCoeffs()
{

    if (this->updated())
    {
        return;
    }

    const label patchI = patch().index();

    const compressible::turbulenceModel& turbModel =
        db().lookupObject<compressible::turbulenceModel>
        (
            IOobject::groupName
            (
                compressible::turbulenceModel::propertiesName,
                internalField().group()
            )
        );

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    const scalarField alphap(turbModel.alphaEff(patchI));


    refValue() = Tinf_;
    refGrad() = 0.0;
    
    valueFraction() =
        1.0/(1.0 + alphap*patch().deltaCoeffs()*patch().magSf()/max(mag(phip), SMALL));

    mixedFvPatchField<scalar>::updateCoeffs();
}


void Foam::fixedEnthalpyFluxTemperatureFvPatchScalarField::
write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeEntry(os, "phi", phiName_);
    writeEntry(os, "rho", rhoName_);
    writeEntry(os, "Tinf", Tinf_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedEnthalpyFluxTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //
