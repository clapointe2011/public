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

#include "injectionModelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

injectionModelList::injectionModelList(surfaceFilmRegionModel& film)
:
    PtrList<injectionModel>(),
    filmSubModelBase(film)
{}


injectionModelList::injectionModelList
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    PtrList<injectionModel>(),
    filmSubModelBase
    (
        "injectionModelList",
        film,
        dict,
        "injectionModelList",
        "injectionModelList"
    ),
    massInjected_(film.intCoupledPatchIDs().size(), 0.0)
{
    const wordList activeModels(dict.lookup("injectionModels"));

    wordHashSet models;
    forAll(activeModels, i)
    {
        models.insert(activeModels[i]);
    }

    Info<< "    Selecting film injection models" << endl;
    if (models.size() > 0)
    {
        this->setSize(models.size());

        label i = 0;
        forAllConstIter(wordHashSet, models, iter)
        {
            const word& model = iter.key();
            set(i, injectionModel::New(film, dict, model));
            i++;
        }
    }
    else
    {
        Info<< "        none" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

injectionModelList::~injectionModelList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void injectionModelList::correct
(
    scalarField& availableMass,
    volScalarField& massToInject,
    volScalarField& diameterToInject
)
{
    // Correct models that accumulate mass and diameter transfers
    forAll(*this, i)
    {
        injectionModel& im = operator[](i);
        im.correct(availableMass, massToInject, diameterToInject);
    }

    // Push values to boundaries ready for transfer to the primary region
    massToInject.correctBoundaryConditions();
    diameterToInject.correctBoundaryConditions();

    const labelList& patchIDs = film().intCoupledPatchIDs();

    forAll(patchIDs, i)
    {
        label patchi = patchIDs[i];
        massInjected_[i] =
            massInjected_[i] + sum(massToInject.boundaryField()[patchi]);
    }
}


void injectionModelList::info(Ostream& os)
{
    const polyBoundaryMesh& pbm = film().regionMesh().boundaryMesh();

    scalar injectedMass = 0;
    scalarField patchInjectedMasses
    (
        pbm.size() - film().regionMesh().globalData().processorPatches().size(),
        0
    );

    forAll(*this, i)
    {
        const injectionModel& im = operator[](i);
        injectedMass += im.injectedMassTotal();
        im.patchInjectedMassTotals(patchInjectedMasses);
    }

    os  << indent << "injected mass      = " << injectedMass << nl;

    forAll(patchInjectedMasses, patchi)
    {
        if (mag(patchInjectedMasses[patchi]) > vSmall)
        {
            os  << indent << indent << "from patch " << pbm[patchi].name()
                << " = " << patchInjectedMasses[patchi] << nl;
        }
    }

    scalarField mass0(massInjected_.size(), 0);
    this->getBaseProperty("massInjected", mass0);

    scalarField mass(massInjected_);
    Pstream::listCombineGather(mass, plusEqOp<scalar>());
    mass += mass0;

    const labelList& patchIDs = film().intCoupledPatchIDs();

    forAll(patchIDs, i)
    {
        label patchi = patchIDs[i];
        Info<< indent << "  - patch: " << pbm[patchi].name() << ": "
            << mass[i] << endl;
    }

    if (film().time().writeTime())
    {
        setBaseProperty("massInjected", mass);
        massInjected_ = 0.0;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
