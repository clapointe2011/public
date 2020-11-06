/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "FLAMES_buoyantWHALE.H"
#include "fvOptions.H"
#include "uniformDimensionedFields.H"
#include "zeroGradientFvPatchFields.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void FLAMES_buoyantWHALE<BasicTurbulenceModel>::correctNut()
{
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);
    BasicTurbulenceModel::correctNut();
}

template<class BasicTurbulenceModel>
tmp<volScalarField> FLAMES_buoyantWHALE<BasicTurbulenceModel>::Gcoef() const
{
    const uniformDimensionedVectorField& g =
        this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");

    return
        Cg_*this->alpha_/(this->k_ + dimensionedScalar(this->k_.dimensions(),VSMALL))
	   *(g & fvc::grad(this->rho_))*this->nut_/this->Prt_;
}

template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> FLAMES_buoyantWHALE<BasicTurbulenceModel>::kSource() const
{

    const uniformDimensionedVectorField& g =
        this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");

    if (mag(g.value()) > small)
    {
        return -fvm::SuSp(Gcoef(), this->k_);
    }
    else
    {
        return tmp<fvScalarMatrix>
        (
            new fvScalarMatrix
            (
                k_,
                dimVolume*this->rho_.dimensions()*k_.dimensions()
                /dimTime
            )
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
FLAMES_buoyantWHALE<BasicTurbulenceModel>::FLAMES_buoyantWHALE
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    LESeddyViscosity<BasicTurbulenceModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    Cw_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw",
            this->coeffDict_,
            0.325
        )
    ),

    Cg_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cg",
            this->coeffDict_,
            1.0
        )
    )
{
    bound(k_, this->kMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool FLAMES_buoyantWHALE<BasicTurbulenceModel>::read()
{
    if (LESeddyViscosity<BasicTurbulenceModel>::read())
    {
        Cw_.readIfPresent(this->coeffDict());
        Cg_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicTurbulenceModel>
tmp<volScalarField> FLAMES_buoyantWHALE<BasicTurbulenceModel>::epsilon() const
{
    return volScalarField::New
    (
        IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
        this->Ce_*k()*sqrt(k())/this->delta()
    );
}

template<class BasicTurbulenceModel>
void FLAMES_buoyantWHALE<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    tmp<volTensorField> tgradU(fvc::grad(U));
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    volSymmTensorField Sij = symm(tgradU());
    volScalarField SuSu = Sij && Sij;
    volTensorField gij = tgradU() & tgradU();
    volTensorField Sd = dev(gij) - skew(gij);
    volScalarField SdSd = Sd && Sd;

    tmp<volScalarField> tS
    (
        new volScalarField
        (
            IOobject
            (
                "tdelta",
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            this->mesh_,
            dimensionedScalar("VSMALL",dimLength, VSMALL),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    volScalarField& S = tS.ref();
    S.ref() = this->delta();
    S.ref() -= min(S.ref());
    S.ref() /= max(S.ref());
    S.correctBoundaryConditions();

    nut = S*sqr(this->Cw_*this->delta())*SdSd*sqrt(SdSd)/(sqr(SuSu)*sqrt(SuSu)
         + SdSd*sqrt(sqrt(SdSd)) + dimensionedScalar("SMALL",dimTime/pow3(sqr(dimTime)),SMALL));
 
    correctNut();

    LESeddyViscosity<BasicTurbulenceModel>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));
    volScalarField G(this->GName(), nut*(tgradU() && dev(twoSymm(tgradU()))));
    tgradU.clear();

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha*rho*G
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, k_)
      - fvm::Sp(this->Ce_*alpha*rho*sqrt(k_)/this->delta(), k_)
      + kSource()
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
