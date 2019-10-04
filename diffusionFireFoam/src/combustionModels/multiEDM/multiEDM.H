/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
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
    Foam::combustionModels::multiEDM
Description
    Eddy dissipation model based on the principle of mixed is burnt.
    Adapted based on https://biblio.ugent.be/publication/8601094 by CL
SourceFiles
    multiEDM.C
\*---------------------------------------------------------------------------*/

#ifndef multiEDM_H
#define multiEDM_H
#include "singleStepCombustion.H"

#include "LESModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                  Class multiEDM Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo, class ThermoType>
class multiEDM
:
    public singleStepCombustion <ReactionThermo, ThermoType>
{
    // Private data

    // Private Member Functions

        //- Return the reciprocal of the turbulent mixing timescale
        tmp<volScalarField> rtAdv() const;

        //- Return the reciprocal of the turbulent diffusion timescale
        tmp<volScalarField> rtDiff() const;
        
        //- Return the reciprocal of the buoyant acceleration timescale
        tmp<volScalarField> rtBuoyAcc() const;
        
        //- Return the reciprocal of the chemistry timescale
        tmp<volScalarField> rtChem() const;
        
        //- Return the reciprocal of the flame timescale
        tmp<volScalarField> rtFlame() const;

        //- Disallow copy construct
        multiEDM(const multiEDM&);

        //- Disallow default bitwise assignment
        void operator=(const multiEDM&);

       
public:

    //- Runtime type information
    TypeName("multiEDM");


    // Constructors

        //- Construct from components
        multiEDM
        (
            const word& modelType,
            ReactionThermo& thermo,
            const compressibleTurbulenceModel& turb,
            const word& combustionProperties
        );


    //- Destructor
    virtual ~multiEDM();


    // Member Functions

        // Evolution

            //- Correct combustion rate
            virtual void correct();


        // I-O

            //- Update properties
            virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "multiEDM.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
