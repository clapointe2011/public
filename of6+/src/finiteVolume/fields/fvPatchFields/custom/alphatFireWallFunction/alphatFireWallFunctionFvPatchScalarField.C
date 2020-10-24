/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

#include "alphatFireWallFunctionFvPatchScalarField.H"
#include "LESModel.H"
#include "basicThermo.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"
//#include "muSgsBuoyantWallFunctionFvPatchScalarField.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

/*
scalar alphatFireWallFunctionFvPatchScalarField::maxExp_ = 50.0;
scalar alphatFireWallFunctionFvPatchScalarField::tolerance_ = 0.01;
label alphatFireWallFunctionFvPatchScalarField::maxIters_ = 10;
*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

/*
void alphatFireWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorIn
        (
            "alphatFireWallFunctionFvPatchScalarField::checkType()"
        )
            << "Patch type for patch " << patch().name() << " must be wall\n"
            << "Current patch type is " << patch().type() << nl
            << exit(FatalError);
    }
}
*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alphatFireWallFunctionFvPatchScalarField::
alphatFireWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    QcFlame_(18000.0),
    QcThreshold_(2000.0)
{
    //checkType();
    //read();
}


alphatFireWallFunctionFvPatchScalarField::
alphatFireWallFunctionFvPatchScalarField
(
    const alphatFireWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    QcFlame_(ptf.QcFlame_),
    QcThreshold_(ptf.QcThreshold_)
{}


alphatFireWallFunctionFvPatchScalarField::
alphatFireWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    QcFlame_(dict.lookupOrDefault<scalar>("QcFlame", 18000.0)),
    QcThreshold_(dict.lookupOrDefault<scalar>("QcThreshold", 2000.0))
{
//    checkType();
    //read();
}


alphatFireWallFunctionFvPatchScalarField::
alphatFireWallFunctionFvPatchScalarField
(
    const alphatFireWallFunctionFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    QcFlame_(tppsf.QcFlame_),
    QcThreshold_(tppsf.QcThreshold_)
{
//    checkType();
}


alphatFireWallFunctionFvPatchScalarField::
alphatFireWallFunctionFvPatchScalarField
(
    const alphatFireWallFunctionFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    QcFlame_(tppsf.QcFlame_),
    QcThreshold_(tppsf.QcThreshold_)
{
//    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
/*
void alphatFireWallFunctionFvPatchScalarField::read()
{

    const LESModel& sgs = db().lookupObject<LESModel>("LESProperties");
    beta_ = readScalar(sgs.thermo().lookup("beta"));
    const IOdictionary& environmentalProperties =
        db().lookupObject<IOdictionary>
        (
            "environmentalProperties"
        );
    dimensionedVector g(environmentalProperties.lookup("g"));
    magG_ = mag(g).value();
}
*/

void alphatFireWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{

    // Get info from the SGS model
    const compressible::turbulenceModel& turbModel =
        db().lookupObject<compressible::turbulenceModel>
        (
            IOobject::groupName
            (
                compressible::turbulenceModel::propertiesName,
                internalField().group()
            )
        );

    const basicThermo& thermo = db().lookupObject<basicThermo>
    (
        "thermophysicalProperties"
    );

    // Field data
    const label patchI = patch().index();

    const scalarField& alphaw = turbModel.alpha()().boundaryField()[patchI];
    scalarField& alphatw = *this;

//    const scalarField& rhow = sgs.rho().boundaryField()[patchI];

    const fvPatchScalarField& Tw = thermo.T().boundaryField()[patchI];

    const scalarField Cpw = thermo.Cp()().boundaryField()[patchI];
//    const scalarField Cpw = thermo.Cp(Tw, patchI);

    const scalarField T(Tw.patchInternalField());

    const scalarField magGradTw(max(mag(Tw.snGrad()), VSMALL));
//    const scalarField magGradTw = mag(Tw.snGrad());

    const scalarField gradTw(Tw.snGrad());


//    const scalarField& ry = patch().deltaCoeffs();

    //const scalarField& phiw = db().lookupObject<surfaceScalarField>("phi").boundaryField()[patchI];
    const scalarField& phiw = patch().lookupPatchField<surfaceScalarField, scalar>("phi");


    // scalar maxFuelMassFlux = 0.0;
    // Populate boundary values
    forAll(alphatw, faceI)
    {
        // Convective heat flux based on laminar thermal diffusivity (use to identify flame location)
        //scalar qwL =  alphaw[faceI]*Cpw[faceI]*magGradTw[faceI];
        scalar qwL =  - alphaw[faceI]*Cpw[faceI]*gradTw[faceI];  //positive when heating the wall
        scalar qw = 0.0;
        scalar alphaEff = alphaw[faceI];
        scalar fuelMassFlux = - phiw[faceI]/patch().magSf()[faceI]*2.5*1000.0; //convert to g/m2/s 
        // maxFuelMassFlux = max(maxFuelMassFlux,fuelMassFlux);

        if (qwL <= VSMALL) //Tw > Tg
        {
            alphatw[faceI] = 0.0;
        }
        else 
        {
            if (fuelMassFlux < 0.1)
            {
                qw = min(max(0,qwL),QcThreshold_)/QcThreshold_*QcFlame_;
            }
            else
            {
                scalar exponent = fuelMassFlux/scalar(10.0);
                exponent = min(50.0,exponent);
                qw = QcFlame_ * (fuelMassFlux/10.0/(Foam::exp(exponent)-scalar(1)));
            }
    
            alphaEff = qw/Cpw[faceI]/max(SMALL,-gradTw[faceI]);
            //alphatw[faceI] = min(1.0, max(0.0, alphaEff - alphaw[faceI]));     
            alphatw[faceI] = min(1.0, alphaEff - alphaw[faceI]);     
        }
 
        if (debug)
        {
            Info<< "    alphaEff       = " << alphaEff << nl
                << "    alphaw         = " << alphaw[faceI] << nl
                << "    alphatw      = " << alphatw[faceI] << nl
                << "    Tw             = " << Tw[faceI] << nl
                << "    T              = "  << T[faceI] << nl
                << "    magGradTw      = "  << magGradTw[faceI] << nl
                << "    Cpw            = "  << Cpw[faceI] << nl
                << "    qwL            = "  << qwL << nl
                << "    qw             = "  << qw << nl
                << "    phi            = "  << phiw[faceI] << nl
                << "    fuelMassFlux[g/m2/s]     = "  << fuelMassFlux << nl
                << endl;
        }

    }
    // reduce(maxFuelMassFlux, maxOp<scalar>());
    // Info << "alphatFireWallFunctionFvPatchScalarField::maxFuelMassFlux[g/m2/s]: " << tab << db().time().timeName() << tab << maxFuelMassFlux << endl;

/*
//  Grab the muSgs patch field using generic/base type
    const fvPatchScalarField& muSgsPatchField =
       patch().lookupPatchField<volScalarField, scalar>("muSgs"); 

    // Perform the type checking
    if (!isA<muSgsBuoyantWallFunctionFvPatchScalarField>(muSgsPatchField))
    {
        FatalErrorIn("alphatFireWallFunctionFvPatchScalarField::evaluate()")
            << "Invalid boundary condition for muSgs" << nl
            << "use muSgsBuoyantWallFunction" << nl
            << endl << abort(FatalError);
    }

//    muSgsBuoyantWallFunctionFvPatchScalarField& muSgsw = 
//       const_cast<muSgsBuoyantWallFunctionFvPatchScalarField&>(muSgsPatchField); 

    const muSgsBuoyantWallFunctionFvPatchScalarField& muSgsw = 
       refCast<const muSgsBuoyantWallFunctionFvPatchScalarField>(muSgsPatchField); 

    muSgsBuoyantWallFunctionFvPatchScalarField& muSgs = 
       const_cast<muSgsBuoyantWallFunctionFvPatchScalarField&>(muSgsw); 

    muSgs.evaluateInAlphaSgs();
*/
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    alphatFireWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
