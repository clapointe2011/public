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
Class
    Foam::dynamicRefineMotionSolverFvMesh
Description
    A fvMesh with built-in refinement.
    Determines which cells to refine/unrefine and does all in update().
        // How often to refine
        refineInterval  1;
        // Field to be refinement on
        field           alpha.water;
        // Refine field in between lower..upper
        lowerRefineLevel 0.001;
        upperRefineLevel 0.999;
        // If value < unrefineLevel (default=great) unrefine
        // unrefineLevel   10;
        // Have slower than 2:1 refinement
        nBufferLayers   1;
        // Refine cells only up to maxRefinement levels
        maxRefinement   2;
        // Stop refinement if maxCells reached
        maxCells        200000;
        // Flux field and corresponding velocity field. Fluxes on changed
        // faces get recalculated by interpolating the velocity. Use 'none'
        // on surfaceScalarFields that do not need to be reinterpolated, use
        // NaN to detect use of mapped variable
        correctFluxes
        (
            (phi none)  // NaN)   // none)
            (nHatf none)   // none)
            (rho*phi none)   // none)
            (ghf none)  // NaN)   // none)
        );
        // Write the refinement level as a volScalarField
        dumpLevel       true;
SourceFiles
    dynamicRefineMotionSolverFvMesh.C

Authoers

    Incorporates hexRef wrapper for multi hexRefs from
    https://github.com/ElsevierSoftwareX/SOFTX_2018_143

    Reference:
    Rettenmaier, D., Deising, D., Ouedraogo, Y., Gjonaj, E., De Gersem, H., Bothe, D., Tropea, C., Marschall, H.
    "Load balanced 2D and 3D adaptive mesh refinement in OpenFOAM" (2019),
    SoftwareX, volume 10, 2352-7110,

    Multifield AMR and mesh motion written by CL

\*---------------------------------------------------------------------------*/

#ifndef dynamicRefineMotionSolverFvMesh_H
#define dynamicRefineMotionSolverFvMesh_H

#include "dynamicFvMesh.H"
#include "hexRef.H" //CL use hexRef from https://github.com/ElsevierSoftwareX/SOFTX_2018_143
#include "PackedBoolList.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class motionSolver;

/*---------------------------------------------------------------------------*\
                     Class dynamicRefineMotionSolverFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicRefineMotionSolverFvMesh
:
    public dynamicFvMesh
{
protected:

        //- Mesh cutting engine
        autoPtr<hexRef> meshCutter_;

        //- Dump cellLevel for postprocessing
        Switch dumpLevel_;

        //- Fluxes to map
        HashTable<word> correctFluxes_;

        //- Number of refinement/unrefinement steps done so far.
        label nRefinementIterations_;

        //- Protected cells (usually since not hexes)
        PackedBoolList protectedCell_;

    // Protected Member Functions

        //CL for mesh motion
        autoPtr<motionSolver> motionPtr_;
        
        //- Optional list of vectorFields to update for mesh motion
        //  For modern solvers using Uf and correctPhi to update the flux
        //  after motion it is not necessary to specify a "velocityFields" list
        velocityMotionCorrection velocityMotionCorrection_;

        //CL for mesh refinement

        //- Count set/unset elements in packedlist.
        static label count(const PackedBoolList&, const unsigned int);

        //- Calculate cells that cannot be refined since would trigger
        //  refinement of protectedCell_ (since 2:1 refinement cascade)
        void calculateProtectedCells(PackedBoolList& unrefineableCell) const;

        //- Read the projection parameters from dictionary
        void readDict();

        //- Refine cells. Update mesh and fields.
        autoPtr<mapPolyMesh> refine(const labelList&);

        //- Unrefine cells. Gets passed in centre points of cells to combine.
        autoPtr<mapPolyMesh> unrefine(const labelList&);

        // Selection of cells to un/refine

            //- Calculates approximate value for refinement level so
            //  we don't go above maxCell
            scalar getRefineLevel
            (
                const label maxCells,
                const label maxRefinement,
                const scalar refineLevel,
                const scalarField&
            ) const;

            //- Get per cell max of connected point
            scalarField maxPointField(const scalarField&) const;

            //- Get point max of connected cell
            scalarField maxCellField(const volScalarField&) const;

            scalarField cellToPoint(const scalarField& vFld) const;

            scalarField error
            (
                const scalarField& fld,
                const scalar minLevel,
                const scalar maxLevel
            ) const;

            //- Select candidate cells for refinement
            virtual void selectRefineCandidates
            (
                const scalar lowerRefineLevel,
                const scalar upperRefineLevel,
                const scalarField& vFld,
                PackedBoolList& candidateCell
            ) const;

            //- Subset candidate cells for refinement
            virtual labelList selectRefineCells
            (
                const label maxCells,
                const label maxRefinement,
                const PackedBoolList& candidateCell
            ) const;

            //- Select points that can be unrefined.
            virtual labelList selectUnrefinePoints
            (
                const scalar unrefineLevel,
                const PackedBoolList& markedCell,
                const scalarField& pFld
            ) const;

            //- Extend markedCell with cell-face-cell.
            void extendMarkedCells(PackedBoolList& markedCell) const;

            //- Check all cells have 8 anchor points
            void checkEightAnchorPoints
            (
                PackedBoolList& protectedCell,
                label& nProtected
            ) const;

private:

        //- Disallow default bitwise copy construct
        dynamicRefineMotionSolverFvMesh(const dynamicRefineMotionSolverFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicRefineMotionSolverFvMesh&);

public:

    //- Runtime type information
    TypeName("dynamicRefineMotionSolverFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit dynamicRefineMotionSolverFvMesh(const IOobject& io);


    //- Destructor
    virtual ~dynamicRefineMotionSolverFvMesh();


    // Member Functions

        //- Direct access to the refinement engine
        const hexRef& meshCutter() const
        {
            return meshCutter_;
        }

        //- Cells which should not be refined/unrefined
        const PackedBoolList& protectedCell() const
        {
            return protectedCell_;
        }

        //- Cells which should not be refined/unrefined
        PackedBoolList& protectedCell()
        {
            return protectedCell_;
        }
        
        //- Return the motionSolver -- CL
        const motionSolver& motion() const;

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Dummy balance function
        virtual bool redistribute();

        //- Dummy balance function
        virtual bool balance();

    // Writing

        //- Write using given format, version and compression
        virtual bool writeObject
        (
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp,
            const bool valid
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

