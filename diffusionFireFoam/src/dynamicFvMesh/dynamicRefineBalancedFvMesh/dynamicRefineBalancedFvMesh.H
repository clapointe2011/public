/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 Tyler Voskuilen
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
    Foam::dynamicRefineBalancedFvMesh

Description
    A fvMesh with built-in refinement and run-time load balancing.

SourceFiles
    dynamicRefineBalancedFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicRefineBalancedFvMesh_H
#define dynamicRefineBalancedFvMesh_H

#include "dynamicRefineFvMesh.H"
#include "decompositionMethod.H"
#include "fvMeshDistribute.H"
#include "mapDistributePolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class dynamicRefineBalancedFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicRefineBalancedFvMesh
:
    public dynamicRefineFvMesh
{
private:

        //- Disallow default bitwise copy construct
        dynamicRefineBalancedFvMesh(const dynamicRefineBalancedFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicRefineBalancedFvMesh&);

        label topParentID(label p);
public:

    //- Runtime type information
    TypeName("dynamicRefineBalancedFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit dynamicRefineBalancedFvMesh(const IOobject& io);


    //- Destructor
    virtual ~dynamicRefineBalancedFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Template to update all volField boundaries
        template<class Type> void correctBoundaries();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "dynamicRefineBalancedFvMeshTemplates.C"
#endif

#endif

// ************************************************************************* //
