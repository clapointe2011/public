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

#include "fixedIncidentRadiationFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "constants.H"
#include "radiationModel.H"
#include "fvDOM.H"
#include "absorptionEmissionModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedIncidentRadiationFvPatchScalarField::
fixedIncidentRadiationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), "undefined", "undefined", "undefined-K"),
    QrIncident_(p.size(), 0.0)
{}


fixedIncidentRadiationFvPatchScalarField::
fixedIncidentRadiationFvPatchScalarField
(
    const fixedIncidentRadiationFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(psf, p, iF, mapper),
    temperatureCoupledBase(patch(), psf),
    QrIncident_(psf.QrIncident_)
{}


fixedIncidentRadiationFvPatchScalarField::
fixedIncidentRadiationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    QrIncident_("QrIncident", dict, p.size())
{
    if (dict.found("value") && dict.found("gradient"))
    {
        fvPatchField<scalar>::operator=(Field<scalar>("value", dict, p.size()));
        gradient() = Field<scalar>("gradient", dict, p.size());
    }
    else
    {
        // Still reading so cannot yet evaluate. Make up a value.
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


fixedIncidentRadiationFvPatchScalarField::
fixedIncidentRadiationFvPatchScalarField
(
    const fixedIncidentRadiationFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(psf, iF),
    temperatureCoupledBase(patch(), psf),
    QrIncident_(psf.QrIncident_)
{}


fixedIncidentRadiationFvPatchScalarField::
fixedIncidentRadiationFvPatchScalarField
(
    const fixedIncidentRadiationFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    temperatureCoupledBase(patch(), ptf),
    QrIncident_(ptf.QrIncident_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void fixedIncidentRadiationFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    m(*this,*this);
}


void fixedIncidentRadiationFvPatchScalarField::rmap
(
    const fvPatchScalarField& psf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(psf, addr);

    const fixedIncidentRadiationFvPatchScalarField& thftpsf =
        refCast<const fixedIncidentRadiationFvPatchScalarField>
        (
            psf
        );

    QrIncident_.rmap(thftpsf.QrIncident_, addr);
}


void fixedIncidentRadiationFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    scalarField intFld(patchInternalField());
    
    const radiationModel& radiation =
        db().lookupObject<radiationModel>("radiationProperties");

    scalarField temissivity
    (
        radiation.absorptionEmission().e()().boundaryField()
        [
            //nbrPatch.index()
            patch().index()
        ]
    );

    gradient() =
        temissivity*
        (
            QrIncident_
          - physicoChemical::sigma.value()*pow4(intFld)  //use internal cell T         
        )/kappa(*this);



    fixedGradientFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        scalar Qr = gSum(kappa(*this)*gradient()*patch().magSf());
        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " -> "
            << " radiativeFlux:" << Qr
            << " walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }
}


void fixedIncidentRadiationFvPatchScalarField::write
(
    Ostream& os
) const
{
    fixedGradientFvPatchScalarField::write(os);
    temperatureCoupledBase::write(os);
    writeEntry(os,"QrIncident",QrIncident_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    fixedIncidentRadiationFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
