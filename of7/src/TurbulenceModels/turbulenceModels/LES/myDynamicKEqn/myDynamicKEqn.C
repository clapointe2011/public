/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "myDynamicKEqn.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::tmp<Foam::volScalarField> myDynamicKEqn<BasicTurbulenceModel>::KK() const
{
    return max
    (
        0.5*(filter_(magSqr(this->U_)) - magSqr(filter_(this->U_))),
        dimensionedScalar(sqr(dimVelocity), small)
    );
}


template<class BasicTurbulenceModel>
volScalarField myDynamicKEqn<BasicTurbulenceModel>::Ck
(
    const volSymmTensorField& D,
    const volScalarField& KK
) const
{
    const volSymmTensorField LL
    (
        dev(filter_(sqr(this->U_)) - (sqr(filter_(this->U_))))
    );

    const volSymmTensorField MM
    (
        this->delta()*(filter_(sqrt(k_)*D) - 2*sqrt(KK + filter_(k_))*filter_(D))
    );

    const volScalarField Ck
    (
       fvc::average(0.5*(LL && MM))
      /(
           fvc::average(magSqr(MM))
         + dimensionedScalar(sqr(MM.dimensions()), vSmall)
       )
    );

    tmp<volScalarField> tfld = Ck;
    return tfld();
}


template<class BasicTurbulenceModel>
volScalarField myDynamicKEqn<BasicTurbulenceModel>::Ce
(
    const volSymmTensorField& D,
    const volScalarField& KK
) const
{

    volScalarField mm =
        pow(KK + filter_(k_), 1.5)/(2*this->delta()) - filter_(pow(k_, 1.5))/this->delta();

    volScalarField ee =
        2*this->delta()*Ck(D, KK)
       *(
           filter_(sqrt(k_)*magSqr(D))
         - 2*sqrt(KK + filter_(k_))*magSqr(filter_(D))
        );

    const volScalarField Ce
    (
        fvc::average(ee*mm)
       /(
            fvc::average(magSqr(mm))
	  + dimensionedScalar(sqr(mm.dimensions()), vSmall)
        )
    );

    tmp<volScalarField> tfld = Ce;
    return tfld();
}


template<class BasicTurbulenceModel>
volScalarField myDynamicKEqn<BasicTurbulenceModel>::Ce() const
{
    const volSymmTensorField D(dev(symm(fvc::grad(this->U_))));
    return Ce(D, KK());
}


template<class BasicTurbulenceModel>
void myDynamicKEqn<BasicTurbulenceModel>::correctNut
(
    const volSymmTensorField& D,
    const volScalarField& KK
)
{
    //CL from https://github.com/AlbertoPa/dynamicSmagorinsky/blob/master/dynamicSmagorinsky.C
    //CL bound nut so that nuEff doesn't go negative
    this->nut_ = max(Ck(D, KK)*sqrt(k_)*this->delta(), -this->nu());
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
void myDynamicKEqn<BasicTurbulenceModel>::correctNut()
{
    correctNut(symm(fvc::grad(this->U_)), KK());
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> myDynamicKEqn<BasicTurbulenceModel>::kSource() const
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
myDynamicKEqn<BasicTurbulenceModel>::myDynamicKEqn
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

    simpleFilter_(this->mesh_),
    filterPtr_(LESfilter::New(this->mesh_, this->coeffDict())),
    filter_(filterPtr_())
{
    bound(k_, this->kMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool myDynamicKEqn<BasicTurbulenceModel>::read()
{
    if (LESeddyViscosity<BasicTurbulenceModel>::read())
    {
        filter_.read(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
tmp<volScalarField> myDynamicKEqn<BasicTurbulenceModel>::epsilon() const
{
    return volScalarField::New
    (
        IOobject::groupName("epsilon", this->alphaRhoPhi_.group()),
        Ce()*k()*sqrt(k())/this->delta()
    );
}


template<class BasicTurbulenceModel>
void myDynamicKEqn<BasicTurbulenceModel>::correct()
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
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    LESeddyViscosity<BasicTurbulenceModel>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volSymmTensorField D(dev(symm(tgradU())));
    const volScalarField G(this->GName(), 2.0*nut*(tgradU() && D));
    tgradU.clear();

    const volScalarField KK(this->KK());

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
    ==
        alpha*rho*G
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, k_)
      - fvm::Sp(Ce(D, KK)*alpha*rho*sqrt(k_)/this->delta(), k_)
      + kSource()
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut(D, KK);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
