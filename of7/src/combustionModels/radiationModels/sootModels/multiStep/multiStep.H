/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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
    Foam::radiationModelsModels::multiStep

Group
    grpRadiationSootSubModels

Description
    Multi-step soot model adapted for foam-7 by Caelan Lapointe.
    Orig. Authhor : Min Zhang in "Implementation of a two-equation soot
    model for sprayFoam," Proceedings of CFD with OpenSource Software, 2019,
    Edited by Nilsson H. http://dx.doi.org/10.17196/OS_CFD#YEAR_2019 

SourceFiles
    multiStep.C

\*---------------------------------------------------------------------------*/

#ifndef multiStep_H
#define multiStep_H

#include "interpolationLookUpTable.H"
#include "sootModel.H"
#include "HashTable.H"
#include "fluidThermo.H"
#include "reactingMixture.H"
//#include "singleStepReactingMixture.H"
#include "turbulentFluidThermoModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiationModels
{
namespace sootModels
{

/*---------------------------------------------------------------------------*\
                    Class multiStep Declaration
\*---------------------------------------------------------------------------*/
template<class ThermoType>
class multiStep
:
    public sootModel
{

    // Static functions

        //- Check mixture in thermo
        static const reactingMixture<ThermoType>& checkThermo
        (
            const fluidThermo&
        );

    // Private data
        volScalarField soot_;

         //- soot model parameters 
        volScalarField PDSoot_;
        volScalarField dNSootdt_;
        volScalarField SootVF_;
        volScalarField SSoot_;
        volScalarField dMSootdt_;
        volScalarField r1_;
        volScalarField r2_;
        volScalarField r3_;
        volScalarField r4_;
        volScalarField r5_;
        volScalarField MSoot_;
        volScalarField NSoot_;
        volScalarField YPrec_;

        //- Soot model dictionary
        dictionary coeffsDict_;

        scalar Sct_;
        scalar MWc_;
        scalar wPrec_;
        scalar rhoS_;
        scalar Ca_;
        scalar Ta_;
        scalar Cb_;
        scalar Cg_;
        scalar Tg_;
        scalar m_;
        scalar q_;
        scalar Cw1_;
        scalar CollEff_;
        scalar Cw2_;
        scalar Tw2_;
        scalar NA_;

        //- Thermo package
        const fluidThermo& thermo_;

public:

    //- Runtime type information
    TypeName("multiStep");

    // Constructors

        //- Construct from components
        multiStep
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& modelType
        );


    //- Destructor
    virtual ~multiStep();


    // Member Functions

        // Edit

            //- Main update/correction routine
            virtual void correct();

        // Access

            //- Return Ysoot
            const volScalarField& soot() const
            {
                return soot_;
            }


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sootModels
} // End namespace radiationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "multiStep.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
