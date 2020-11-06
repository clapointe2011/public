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

#include "buoyantWHALE.H"
#include "fvOptions.H"
#include "uniformDimensionedFields.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void buoyantWHALE<BasicTurbulenceModel>::correctNut()
{
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);
    BasicTurbulenceModel::correctNut();
}

template<class BasicTurbulenceModel>
tmp<volScalarField> buoyantWHALE<BasicTurbulenceModel>::Gcoef() const
{
    const uniformDimensionedVectorField& g =
        this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");

    return
        Cg_*this->alpha_/(this->k_ + dimensionedScalar(this->k_.dimensions(),VSMALL))
	   *(g & fvc::grad(this->rho_))*this->nut_/this->Prt_;
}

template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> buoyantWHALE<BasicTurbulenceModel>::kSource() const
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
buoyantWHALE<BasicTurbulenceModel>::buoyantWHALE
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
bool buoyantWHALE<BasicTurbulenceModel>::read()
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
tmp<volScalarField> buoyantWHALE<BasicTurbulenceModel>::epsilon() const
{
    return volScalarField::New
    (
        IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
        this->Ce_*k()*sqrt(k())/this->delta()
    );
}

template<class BasicTurbulenceModel>
void buoyantWHALE<BasicTurbulenceModel>::correct()
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

    volScalarField magSqrSd = magSqr(dev(symm(tgradU() & tgradU())));
    volScalarField magSqrSu = magSqr(symm(tgradU()));

    nut = sqr(Cw_*this->delta())*
        (
            pow(magSqrSd, 3.0/2.0)
           /(
               (
                   pow(magSqrSu, 5.0/2.0)
                 + pow(magSqrSd, 5.0/4.0)
               )
             + dimensionedScalar
               (
                   "small",
                   dimensionSet(0, 0, -5, 0, 0),
                   small
               )
           )
        );

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
