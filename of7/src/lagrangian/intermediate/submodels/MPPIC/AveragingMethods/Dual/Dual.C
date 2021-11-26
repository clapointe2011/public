/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

#include "Dual.H"
#include "coupledPointPatchField.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::labelList> Foam::AveragingMethods::Dual<Type>::size
(
    const fvMesh& mesh
)
{
    autoPtr<labelList> s(new labelList(2));
    s()[0] = mesh.nCells();
    s()[1] = mesh.nPoints();
    return s;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::AveragingMethods::Dual<Type>::Dual
(
    const IOobject& io,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    AveragingMethod<Type>(io, dict, mesh, size(mesh)),
    volumeCell_(mesh.V()),
    volumeDual_(mesh.nPoints(), 0.0),
    dataCell_(FieldField<Field, Type>::operator[](0)),
    dataDual_(FieldField<Field, Type>::operator[](1))
{
    forAll(this->mesh_.C(), celli)
    {
        List<tetIndices> cellTets =
            polyMeshTetDecomposition::cellTetIndices(this->mesh_, celli);
        forAll(cellTets, tetI)
        {
            const tetIndices& tetIs = cellTets[tetI];
            const triFace triIs = tetIs.faceTriIs(this->mesh_);
            const scalar v = tetIs.tet(this->mesh_).mag();
            volumeDual_[triIs[0]] += v;
            volumeDual_[triIs[1]] += v;
            volumeDual_[triIs[2]] += v;
        }
    }

    mesh.globalData().syncPointData
    (
        volumeDual_,
        plusEqOp<scalar>(),
        mapDistribute::transform()
    );
}


template<class Type>
Foam::AveragingMethods::Dual<Type>::Dual
(
    const Dual<Type>& am
)
:
    AveragingMethod<Type>(am),
    volumeCell_(am.volumeCell_),
    volumeDual_(am.volumeDual_),
    dataCell_(FieldField<Field, Type>::operator[](0)),
    dataDual_(FieldField<Field, Type>::operator[](1))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::AveragingMethods::Dual<Type>::~Dual()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::AveragingMethods::Dual<Type>::syncDualData()
{
    this->mesh_.globalData().syncPointData
    (
        dataDual_,
        plusEqOp<Type>(),
        mapDistribute::transform()
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::AveragingMethods::Dual<Type>::add
(
    const barycentric& coordinates,
    const tetIndices& tetIs,
    const Type& value
)
{
    const triFace triIs(tetIs.faceTriIs(this->mesh_));

    dataCell_[tetIs.cell()] +=
        coordinates[0]*value
      / (0.25*volumeCell_[tetIs.cell()]);

    for(label i = 0; i < 3; i ++)
    {
        dataDual_[triIs[i]] +=
            coordinates[i+1]*value
          / (0.25*volumeDual_[triIs[i]]);
    }
}


template<class Type>
Type Foam::AveragingMethods::Dual<Type>::interpolate
(
    const barycentric& coordinates,
    const tetIndices& tetIs
) const
{
    const triFace triIs(tetIs.faceTriIs(this->mesh_));

    return
        coordinates[0]*dataCell_[tetIs.cell()]
      + coordinates[1]*dataDual_[triIs[0]]
      + coordinates[2]*dataDual_[triIs[1]]
      + coordinates[3]*dataDual_[triIs[2]];
}


template<class Type>
typename Foam::AveragingMethods::Dual<Type>::TypeGrad
Foam::AveragingMethods::Dual<Type>::interpolateGrad
(
    const barycentric& coordinates,
    const tetIndices& tetIs
) const
{
    const triFace triIs(tetIs.faceTriIs(this->mesh_));

    const label celli(tetIs.cell());

    const tensor T
    (
        inv
        (
            tensor
            (
                this->mesh_.points()[triIs[0]] - this->mesh_.C()[celli],
                this->mesh_.points()[triIs[1]] - this->mesh_.C()[celli],
                this->mesh_.points()[triIs[2]] - this->mesh_.C()[celli]
            )
        )
    );

    const vector t( - T.T().x() - T.T().y() - T.T().z());

    const TypeGrad S
    (
        dataDual_[triIs[0]],
        dataDual_[triIs[1]],
        dataDual_[triIs[2]]
    );

    const Type s(dataCell_[celli]);

    return (T & S) + (t*s);
}


template<class Type>
void Foam::AveragingMethods::Dual<Type>::average()
{
    syncDualData();

    AveragingMethod<Type>::average();
}


template<class Type>
void Foam::AveragingMethods::Dual<Type>::average
(
    const AveragingMethod<scalar>& weight
)
{
    syncDualData();

    AveragingMethod<Type>::average(weight);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::AveragingMethods::Dual<Type>::primitiveField() const
{
    return tmp<Field<Type>>(dataCell_);
}


// ************************************************************************* //
