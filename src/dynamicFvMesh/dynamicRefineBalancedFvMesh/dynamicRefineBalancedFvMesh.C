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


\*---------------------------------------------------------------------------*/

#include "dynamicRefineBalancedFvMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicRefineBalancedFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, dynamicRefineBalancedFvMesh, IOobject);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::dynamicRefineBalancedFvMesh::topParentID(label p)
{
    label nextP = meshCutter().history().splitCells()[p].parent_;
    if( nextP < 0 )
    {
        return p;
    }
    else
    {
        return topParentID(nextP);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicRefineBalancedFvMesh::dynamicRefineBalancedFvMesh
(
    const IOobject& io
)
:
    dynamicRefineFvMesh(io)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicRefineBalancedFvMesh::~dynamicRefineBalancedFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicRefineBalancedFvMesh::update()
{
//Part 1 - Call normal update from dynamicRefineFvMesh
    bool hasChanged = dynamicRefineFvMesh::update();

// Part 2 - Load Balancing
    dictionary refineDict
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).subDict("dynamicRefineFvMeshCoeffs")
    );

    Switch enableBalancing = refineDict.lookup("enableBalancing");

    if ( Pstream::parRun() && hasChanged )
    {
        const scalar allowableImbalance =
            readScalar(refineDict.lookup("allowableImbalance"));

        //First determine current level of imbalance - do this for all
        // parallel runs with a changing mesh, even if balancing is disabled
        label nGlobalCells = globalData().nTotalCells();
        scalar idealNCells = scalar(nGlobalCells)/scalar(Pstream::nProcs());
        scalar localImbalance = mag(scalar(nCells()) - idealNCells);
        Foam::reduce(localImbalance, maxOp<scalar>());
        scalar maxImbalance = localImbalance/idealNCells;

        Info<<"Maximum imbalance = " << 100*maxImbalance << " %" << endl;

        //If imbalanced, construct weighted coarse graph (level 0) with node
        // weights equal to their number of subcells. This partitioning works
        // as long as the number of level 0 cells is several times greater than
        // the number of processors.
        if( maxImbalance > allowableImbalance && enableBalancing)
        {
            Info<< "Re-balancing dynamically refined mesh" << endl;

            const labelIOList& cellLevel = meshCutter().cellLevel();
            Map<label> coarseIDmap(100);
            labelList uniqueIndex(nCells(),0);

            label nCoarse = 0;

            forAll(cells(), cellI)
            {
                if( cellLevel[cellI] > 0 )
                {
                    uniqueIndex[cellI] = nCells() + topParentID
                    (
                        meshCutter().history().parentIndex(cellI)
                    );
                }
                else
                {
                    uniqueIndex[cellI] = cellI;
                }

                if( coarseIDmap.insert(uniqueIndex[cellI], nCoarse) )
                {
                    ++nCoarse;
                }
            }

            // Convert to local sequential indexing and calculate coarse
            // points and weights
            labelList localIndex(nCells(),0);
            pointField coarsePoints(nCoarse,vector::zero);
            scalarField coarseWeights(nCoarse,0.0);

            forAll(uniqueIndex, cellI)
            {
                localIndex[cellI] = coarseIDmap[uniqueIndex[cellI]];

                // If 2D refinement (quadtree) is ever implemented, this '3'
                // should be set in general as the number of refinement
                // dimensions.
                label w = (1 << (3*cellLevel[cellI]));

                coarseWeights[localIndex[cellI]] += 1.0;
                coarsePoints[localIndex[cellI]] += C()[cellI]/w;
            }

            //Set up decomposer - a separate dictionary is used here so
            // you can use a simple partitioning for decomposePar and
            // ptscotch for the rebalancing (or any chosen algorithms)
            autoPtr<decompositionMethod> decomposer
            (
                decompositionMethod::New
                (
                    IOdictionary
                    (
                        IOobject
                        (
                            "balanceParDict",
                            time().system(),
                            *this,
                            IOobject::MUST_READ_IF_MODIFIED,
                            IOobject::NO_WRITE
                        )
                    )
                )
            );

            //Pout<<"debug message prior to decomposer().decompose(...)"<<endl;

            labelList finalDecomp = decomposer().decompose
            (
                *this,          //polyMesh
                localIndex,     //labelList
                coarsePoints,   //pointField
                coarseWeights   //scalarField
            );

            //Pout<<"debug mesage after decomposer().decompose(...)"<<endl;

            scalar tolDim = globalMeshData::matchTol_ * bounds().mag();

            fvMeshDistribute distributor(*this, tolDim);

            autoPtr<mapDistributePolyMesh> map =
                distributor.distribute(finalDecomp);

            meshCutter_.distribute(map);

            //Correct values on all cyclic patches
            correctBoundaries<scalar>();
            correctBoundaries<vector>();
            correctBoundaries<sphericalTensor>();
            correctBoundaries<symmTensor>();
            correctBoundaries<tensor>();
        }
    }

    return hasChanged;
}


// ************************************************************************* //
