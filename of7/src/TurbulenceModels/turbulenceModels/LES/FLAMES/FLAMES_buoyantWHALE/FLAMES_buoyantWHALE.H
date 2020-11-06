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

Class
    Foam::LESModels::FLAMES_buoyantWHALE

Description
    One equation eddy-viscosity model

    Eddy viscosity SGS model using a modeled balance equation to simulate the
    behaviour of k.

    Based on WALE model from FM Global fireFoam-dev
SourceFiles
    FLAMES_buoyantWHALE.C

\*---------------------------------------------------------------------------*/

#ifndef FLAMES_buoyantWHALE_H
#define FLAMES_buoyantWHALE_H

#include "LESeddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                           Class FLAMES_buoyantWHALE Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class FLAMES_buoyantWHALE
:
    public LESeddyViscosity<BasicTurbulenceModel>
{
protected:

    // Protected data

        // Fields

            volScalarField k_;


        // Model constants

            dimensionedScalar Cw_;
            dimensionedScalar Cg_;

    // Protected Member Functions
        tmp<volScalarField> Gcoef() const;
        virtual void correctNut();
        virtual tmp<fvScalarMatrix> kSource() const;

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("FLAMES_buoyantWHALE");


    // Constructors

        //- Constructor from components
        FLAMES_buoyantWHALE
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );

        //- Disallow default bitwise copy construction
        FLAMES_buoyantWHALE(const FLAMES_buoyantWHALE&) = delete;


    //- Destructor
    virtual ~FLAMES_buoyantWHALE()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return sub-grid disipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return volScalarField::New
            (
                "DkEff",
                this->nut_ + this->nu()
            );
        }

        //- Correct eddy-Viscosity and related properties
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const FLAMES_buoyantWHALE&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FLAMES_buoyantWHALE.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
