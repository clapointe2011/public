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

#include "filmPyrolysisRadiativeCoupledMixed2FvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "mappedPatchBase.H"
#include "constants.H"
#include "radiationModel.H"
#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const filmPyrolysisRadiativeCoupledMixed2FvPatchScalarField::filmModelType&
filmPyrolysisRadiativeCoupledMixed2FvPatchScalarField::
filmModel() const
{
    HashTable<const filmModelType*> models
        = db().time().lookupClass<filmModelType>();

    forAllConstIter(HashTable<const filmModelType*>, models, iter)
    {
        if (iter()->regionMesh().name() == filmRegionName_)
        {
            return *iter();
        }
    }

    DynamicList<word> modelNames;
    forAllConstIter(HashTable<const filmModelType*>, models, iter)
    {
        modelNames.append(iter()->regionMesh().name());
    }

    FatalErrorInFunction
        << "Unable to locate film region " << filmRegionName_
        << ".  Available regions include: " << modelNames
        << abort(FatalError);

    return **models.begin();
}


const filmPyrolysisRadiativeCoupledMixed2FvPatchScalarField::
pyrolysisModelType&
filmPyrolysisRadiativeCoupledMixed2FvPatchScalarField::
pyrModel() const
{
    HashTable<const pyrolysisModelType*> models =
        db().time().lookupClass<pyrolysisModelType>();

    forAllConstIter(HashTable<const pyrolysisModelType*>, models, iter)
    {
        if (iter()->regionMesh().name() == pyrolysisRegionName_)
        {
            return *iter();
        }
    }

    DynamicList<word> modelNames;
    forAllConstIter(HashTable<const pyrolysisModelType*>, models, iter)
    {
        modelNames.append(iter()->regionMesh().name());
    }


    FatalErrorInFunction
        << "Unable to locate pyrolysis region " << pyrolysisRegionName_
        << ".  Available regions include: " << modelNames
        << abort(FatalError);

    return **models.begin();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

filmPyrolysisRadiativeCoupledMixed2FvPatchScalarField::
filmPyrolysisRadiativeCoupledMixed2FvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), "undefined", "undefined", "undefined-K"),
    radiationCoupledBase(p, "undefined", scalarField::null()),
    filmRegionName_("surfaceFilmProperties"),
    pyrolysisRegionName_("pyrolysisProperties"),
    TnbrName_("undefined-Tnbr"),
    qrNbrName_("undefined-qrNbr"),
    qrName_("undefined-qr"),
    convectiveScaling_(1.0),
    filmDeltaDry_(0.0),
    filmDeltaWet_(0.0),
    emissivity_(p.size(), 0.0)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


filmPyrolysisRadiativeCoupledMixed2FvPatchScalarField::
filmPyrolysisRadiativeCoupledMixed2FvPatchScalarField
(
    const filmPyrolysisRadiativeCoupledMixed2FvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    temperatureCoupledBase(patch(), psf),
    radiationCoupledBase
    (
        p,
        psf.emissivityMethod(),
        psf.emissivity_
    ),
    filmRegionName_(psf.filmRegionName_),
    pyrolysisRegionName_(psf.pyrolysisRegionName_),
    TnbrName_(psf.TnbrName_),
    qrNbrName_(psf.qrNbrName_),
    qrName_(psf.qrName_),
    convectiveScaling_(psf.convectiveScaling_),
    filmDeltaDry_(psf.filmDeltaDry_),
    filmDeltaWet_(psf.filmDeltaWet_),
    emissivity_(psf.emissivity_)
{}


filmPyrolysisRadiativeCoupledMixed2FvPatchScalarField::
filmPyrolysisRadiativeCoupledMixed2FvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    radiationCoupledBase(p, dict),
    filmRegionName_
    (
        dict.lookupOrDefault<word>("filmRegion", "surfaceFilmProperties")
    ),
    pyrolysisRegionName_
    (
        dict.lookupOrDefault<word>("pyrolysisRegion", "pyrolysisProperties")
    ),
    TnbrName_(dict.lookup("Tnbr")),
    qrNbrName_(dict.lookup("qrNbr")),
    qrName_(dict.lookup("qr")),
    convectiveScaling_(dict.lookupOrDefault<scalar>("convectiveScaling", 1.0)),
    filmDeltaDry_
    (
            dict.lookupOrDefault<scalar>("filmDeltaDry", 0.0000)
    ),
    filmDeltaWet_
    (
            dict.lookupOrDefault<scalar>("filmDeltaWet", 0.0002)
    ),
    emissivity_(p.size(), 0.0)
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorInFunction
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}


filmPyrolysisRadiativeCoupledMixed2FvPatchScalarField::
filmPyrolysisRadiativeCoupledMixed2FvPatchScalarField
(
    const filmPyrolysisRadiativeCoupledMixed2FvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    temperatureCoupledBase(patch(), psf),
    radiationCoupledBase
    (
        psf.patch(),
        psf.emissivityMethod(),
        psf.emissivity_
    ),
    filmRegionName_(psf.filmRegionName_),
    pyrolysisRegionName_(psf.pyrolysisRegionName_),
    TnbrName_(psf.TnbrName_),
    qrNbrName_(psf.qrNbrName_),
    qrName_(psf.qrName_),
    convectiveScaling_(psf.convectiveScaling_),
    filmDeltaDry_(psf.filmDeltaDry_),
    filmDeltaWet_(psf.filmDeltaWet_),
    emissivity_(psf.emissivity_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void filmPyrolysisRadiativeCoupledMixed2FvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());

    const label patchi = patch().index();
    const label nbrPatchi = mpp.samplePolyPatch().index();
    const polyMesh& mesh = patch().boundaryMesh().mesh();
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[nbrPatchi];

    scalarField intFld(patchInternalField());

    const filmPyrolysisRadiativeCoupledMixed2FvPatchScalarField&
        nbrField =
        refCast
        <
            const filmPyrolysisRadiativeCoupledMixed2FvPatchScalarField
        >
        (
            nbrPatch.lookupPatchField<volScalarField, scalar>(TnbrName_)
        );

    // Swap to obtain full local values of neighbour internal field
    scalarField nbrIntFld(nbrField.patchInternalField());
    mpp.distribute(nbrIntFld);

    const scalarField K(this->kappa(*this));
    const scalarField nbrK(nbrField.kappa(*this));

    // Swap to obtain full local values of neighbour K*delta
    scalarField KDeltaNbr(nbrK*nbrPatch.deltaCoeffs());
    mpp.distribute(KDeltaNbr);
    scalarField nbrTotalFlux(patch().size(), 0.0); //CL from firefoam-dev
    scalarList Tfilm(patch().size(), 0.0);
    scalarList filmDelta(patch().size(), 0.0);
    scalarList filmConv(patch().size(), 0.0); //CL from firefoam-dev
    scalarList nbrConv(patch().size(), 0.0);
    scalarList alpha(patch().size(), 0.0); //CL from firefoam-dev
    scalarField qr(patch().size(), 0.0);
    scalarList qrCoupled(nbrPatch.size(), 0.0); //CL from firefoam-dev

    const pyrolysisModelType& pyrolysis = pyrModel();
    const filmModelType& film = filmModel();

    label coupledPatchi = -1;
    if (pyrolysisRegionName_ == mesh.name())
    {
        coupledPatchi = patchi;
        if (qrName_ != "none")
        {
            qr = nbrPatch.lookupPatchField<volScalarField, scalar>(qrName_);
            mpp.distribute(qr);
        }
    }
    else if (pyrolysis.primaryMesh().name() == mesh.name())
    {
        coupledPatchi = nbrPatch.index();
        if (qrName_ != "none")
        {
            qr = patch().lookupPatchField<volScalarField, scalar>(qrName_);
        }
    }
    else
    {
        FatalErrorInFunction
            << type() << " condition is intended to be applied to either the "
            << "primary or pyrolysis regions only"
            << exit(FatalError);
    }

    // In solid, qr = qrNbr
    if(qrNbrName_ != "none")
    {
        const label filmPatchi = pyrolysis.nbrCoupledPatchID(film, coupledPatchi);
        const scalarField Qconvw(film.Qconvw(filmPatchi));
        const scalarField nbrQconv(convectiveScaling_*KDeltaNbr*(intFld - nbrIntFld));

        // Obtain film convective heat transfer
	filmConv =
            pyrolysis.mapRegionPatchField
            (
                 film,
                 coupledPatchi,
                 filmPatchi,
                 Qconvw,
                 true
            );

        // Obtain neighbour convective heat transfer
        nbrConv =
            pyrolysis.mapRegionPatchField
            (
                 film,
                 coupledPatchi,
                 filmPatchi,
                 nbrQconv,
                 true
            );

        /*qrCoupled =
            pyrolysis.mapRegionPatchField<scalar>
            (
                film,
                "qin",
                coupledPatchi,
                true
            );
        */

        // Obtain delta
        filmDelta =
            pyrolysis.mapRegionPatchField<scalar>
            (
                film,
                "deltaf",
                coupledPatchi,
                true
            );

        // Obtain alpha
        alpha =
            pyrolysis.mapRegionPatchField<scalar>
            (
                film,
                "alpha",
                coupledPatchi,
                true
            );

        const radiationModel& radiation =
        db().lookupObject<radiationModel>
        (
            "radiationProperties"
        );

        scalarField temissivity
        (
            radiation.absorptionEmission().e()().boundaryField()
            [
                patch().index()
            ]
        );

        forAll(*this,i)
        {
            scalar qConvWeti = -filmConv[i];
            scalar qConvDryi = nbrConv[i];

            scalar qRadWeti = 0.0; // all film absorption takes place in film model
            scalar qRadDryi = qr[i];
            //   -temissivity[i]*qrCoupled[i]
            //   +temissivity[i]*constant::physicoChemical::sigma.value()*pow(operator[](i),4);

            scalar qConvi = alpha[i]*qConvWeti + (1.0-alpha[i])*qConvDryi;
            scalar qRadi  = alpha[i]*qRadWeti  + (1.0-alpha[i])*qRadDryi;

	    nbrTotalFlux[i] = qConvi + qRadi;
            this->refValue()[i] = operator[](i);  // not used
            this->refGrad()[i] = -nbrTotalFlux[i]/K[i];
            this->valueFraction()[i] = 0.0;
        }
    }
    //In fluid, qr = qr
    else
    {
        Tfilm =
            pyrolysis.mapRegionPatchField<scalar>
            (
                film,
                "Tf",
                nbrPatchi,
                true
            );

        filmDelta =
            pyrolysis.mapRegionPatchField<scalar>
            (
                film,
                "deltaf",
                nbrPatchi,
                true
            );

        // do we still need to do mpp.distribute Tfilm
        mpp.distribute(Tfilm);

        mpp.distribute(filmDelta);

        scalarList Twall(patch().size(), 0.0);
        
        // Estimate wetness of the film (1: wet , 0: dry)
        scalarField ratio
        (
           min
           (
               max
               (
                   (filmDelta - filmDeltaDry_)/(filmDeltaWet_ - filmDeltaDry_),
                   scalar(0.0)
               ),
               scalar(1.0)
           )
        );

        forAll(*this, i)
        {
            scalar Twet = min(max(Tfilm[i], 298.15), 378.4);
            scalar Tdry = nbrIntFld[i];

            Twall[i] = ratio[i]*(Twet - Tdry) + Tdry;
        }

        this->refValue() = Twall;
        this->refGrad() = 0.0;   // not used
        this->valueFraction() = 1.0;
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void filmPyrolysisRadiativeCoupledMixed2FvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    writeEntryIfDifferent<word>
    (
        os,
        "filmRegion",
        "surfaceFilmProperties",
        filmRegionName_
    );
    writeEntryIfDifferent<word>
    (
        os,
        "pyrolysisRegion",
        "pyrolysisProperties",
        pyrolysisRegionName_
    );
    writeEntry(os, "Tnbr", TnbrName_);
    writeEntry(os, "qr", qrName_);
    writeEntry(os, "convectiveScaling", convectiveScaling_);
    writeEntry(os, "filmDeltaDry", filmDeltaDry_);
    writeEntry(os, "filmDeltaWet", filmDeltaWet_);
    temperatureCoupledBase::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    filmPyrolysisRadiativeCoupledMixed2FvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
