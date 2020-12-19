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

#include "myReactingOneDim.H"
#include "addToRunTimeSelectionTable.H"
#include "fvm.H"
#include "fvcDiv.H"
#include "fvcVolumeIntegrate.H"
#include "fvcLaplacian.H"
#include "absorptionEmissionModel.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace pyrolysisModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(myReactingOneDim, 0);

addToRunTimeSelectionTable(pyrolysisModel, myReactingOneDim, mesh);
addToRunTimeSelectionTable(pyrolysisModel, myReactingOneDim, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void myReactingOneDim::readmyReactingOneDimControls()
{
    const dictionary& solution = this->solution().subDict("SIMPLE");
    solution.lookup("nNonOrthCorr") >> nNonOrthCorr_;
    time().controlDict().lookup("maxDi") >> maxDiff_;

    coeffs().lookup("minimumDelta") >> minimumDelta_;

    coeffs().lookup("gasHSource") >> gasHSource_;
    coeffs().lookup("qrHSource") >> qrHSource_;
    coeffs().lookup("charOxiSource") >> charOxiSource_;
    useChemistrySolvers_ =
        coeffs().lookupOrDefault<bool>("useChemistrySolvers", true);
}


bool myReactingOneDim::read()
{
    if (pyrolysisModel::read())
    {
        readmyReactingOneDimControls();
        return true;
    }
    else
    {
        return false;
    }
}


bool myReactingOneDim::read(const dictionary& dict)
{
    if (pyrolysisModel::read(dict))
    {
        readmyReactingOneDimControls();
        return true;
    }
    else
    {
        return false;
    }
}

void myReactingOneDim::updateCharOxi()
{
    const label charIndex = solidThermo_->composition().species()["char"];

    scalar rhoChar = solidThermo_->composition().rho(charIndex, 10000, 298);

    //Info << "mw of char is " <<  solidThermo_.composition().W(charIndex);

    // hardcode molecular weight for C and O2
    scalar mWO2 = 32;
    scalar mWChar = 12;
    scalar mWCO2 = 44;

    // hardcode HoC for char [J/kg]
    dimensionedScalar HocChar("HocChar",dimEnergy/dimMass,32.8e6);

    label localPyrolysisFaceI = 0;

    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];
        scalarField& Xcharp = Xchar_.boundaryFieldRef()[patchi];
        scalarField& mCharp = mChar_.boundaryFieldRef()[patchi];
        scalarField& mCharBurntp = mCharBurnt_.boundaryFieldRef()[patchi];
        scalarField& charOxiQdotp = charOxiQdot_.boundaryFieldRef()[patchi];
        scalarField& phiO2p = phiO2_.boundaryFieldRef()[patchi];
        scalarField& phiCO2p = phiCO2_.boundaryFieldRef()[patchi];

        const fvPatch& patch = regionMesh().boundary()[patchi];

        // Get the coupling information from the mappedPatchBase
        const mappedPatchBase& mpp = refCast<const mappedPatchBase>
        (
            patch.patch()
        );
        const polyMesh& nbrMesh = mpp.sampleMesh();
        const fvPatch& nbrPatch = refCast<const fvMesh>
        (
            nbrMesh
        ).boundary()[mpp.samplePolyPatch().index()];

        fvPatchScalarField O2 = nbrPatch.lookupPatchField<volScalarField, scalar>("O2");
        scalarField alpha = nbrPatch.lookupPatchField<volScalarField, scalar>("thermo:alpha");    //for now use molecular alpha
        scalarField O2Int(O2.patchInternalField());
        scalarField alphaDelta(alpha * nbrPatch.deltaCoeffs());

        mpp.distribute(O2);
        mpp.distribute(O2Int);
        mpp.distribute(alphaDelta);

        phiO2p = - alphaDelta * (O2Int - 0.0) * patch.magSf(); //[kg/s] negative
        phiCO2p = - phiO2p * mWCO2 / mWO2;  // positive, same as phiGas

        scalarField deltaMO2(- phiO2p * time_.deltaT().value());   //[kg]

        const scalarField& cellV = regionMesh().V();

        forAll(Xcharp, faceI)
        {
            const labelList& cells = boundaryFaceCells_[localPyrolysisFaceI];
            scalar charVolume = 0.0;
            scalar totalVolume = 0.0;
            forAll(cells, k)
            {
                const label cellI = cells[k];
                scalar X = Ys_[charIndex][cellI] * rho_[cellI] / rhoChar;

                charVolume += X * cellV[cellI];
                totalVolume += cellV[cellI];
            }
            Xcharp[faceI] = charVolume / totalVolume;
            mCharp[faceI] = rhoChar * charVolume;

            scalar charAvai = mCharp[faceI] - mCharBurntp[faceI];
            scalar dmCharBurnt = 0.0;
            if (charAvai < deltaMO2[faceI] / mWO2 * mWChar)
            {
                dmCharBurnt = charAvai;
                phiO2p[faceI] = - dmCharBurnt / mWChar * mWO2;
                phiCO2p[faceI] = dmCharBurnt / mWChar * mWCO2;
            }
            else
            {
                dmCharBurnt = deltaMO2[faceI] / mWO2 * mWChar;
            }
            mCharBurntp[faceI] += dmCharBurnt;

            //compute energy release rate in the first cell [kW/m3]
            charOxiQdot_[cells[0]] = HocChar.value() * dmCharBurnt
                                   / cellV[cells[0]]
                                   / time_.deltaT().value();

            //for HRR diagnostics (sum at the patch is the HRR) [kW]
            charOxiQdotp[faceI] = charOxiQdot_[cells[0]] * cellV[cells[0]];

            localPyrolysisFaceI++;
        }
    }
}


void myReactingOneDim::updateqr()
{
    // Update local qr from coupled qr field
    qr_ == dimensionedScalar("zero", qr_.dimensions(), 0.0);

    // Retrieve field from coupled region using mapped boundary conditions
    qr_.correctBoundaryConditions();

    volScalarField::Boundary& qrBf = qr_.boundaryFieldRef();

    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];

        // qr is positive going in the solid
        // If the surface is emitting the radiative flux is set to zero
        qrBf[patchi] = max(qrBf[patchi], scalar(0));
    }

    const vectorField& cellC = regionMesh().cellCentres();

    tmp<volScalarField> kappa = kappaRad();

    // Propagate qr through 1-D regions
    label localPyrolysisFacei = 0;
    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];

        const scalarField& qrp = qr_.boundaryField()[patchi];
        const vectorField& Cf = regionMesh().Cf().boundaryField()[patchi];

        forAll(qrp, facei)
        {
            const scalar qr0 = qrp[facei];
            point Cf0 = Cf[facei];
            const labelList& cells = boundaryFaceCells_[localPyrolysisFacei++];
            scalar kappaInt = 0.0;
            forAll(cells, k)
            {
                const label celli = cells[k];
                const point& Cf1 = cellC[celli];
                const scalar delta = mag(Cf1 - Cf0);
                kappaInt += kappa()[celli]*delta;
                qr_[celli] = qr0*exp(-kappaInt);
                Cf0 = Cf1;
            }
        }
    }
}


void myReactingOneDim::updatePhiGas()
{
    phiHsGas_ ==  dimensionedScalar("zero", phiHsGas_.dimensions(), 0.0);
    phiGas_ == dimensionedScalar("zero", phiGas_.dimensions(), 0.0);

    const speciesTable& gasTable = solidChemistry_->gasTable();

    forAll(gasTable, gasI)
    {
        tmp<volScalarField> tHsiGas =
            solidChemistry_->gasHs(solidThermo_->p(), solidThermo_->T(), gasI);

        const volScalarField& HsiGas = tHsiGas();

        const volScalarField::Internal& RRiGas =
            solidChemistry_->RRg(gasI);

        surfaceScalarField::Boundary& phiGasBf =
            phiGas_.boundaryFieldRef();

        label totalFaceId = 0;
        forAll(intCoupledPatchIDs_, i)
        {
            const label patchi = intCoupledPatchIDs_[i];

            scalarField& phiGasp = phiGasBf[patchi];
            const scalarField& cellVol = regionMesh().V();

            forAll(phiGasp, facei)
            {
                const labelList& cells = boundaryFaceCells_[totalFaceId++];
                scalar massInt = 0.0;
                forAllReverse(cells, k)
                {
                    const label celli = cells[k];
                    massInt += RRiGas[celli]*cellVol[celli];
                    phiHsGas_[celli] += massInt*HsiGas[celli];
                }

                phiGasp[facei] += massInt;

                if (debug)
                {
                    Info<< " Gas : " << gasTable[gasI]
                        << " on patch : " << patchi
                        << " mass produced at face(local) : "
                        <<  facei
                        << " is : " << massInt
                        << " [kg/s] " << endl;
                }
            }
        }
    }
}


void myReactingOneDim::updateFields()
{
    if (qrHSource_)
    {
        updateqr();
    }

    if (charOxiSource_)
    {
        updateCharOxi();
    }

    updatePhiGas();
}


void myReactingOneDim::updateMesh(const scalarField& mass0)
{
    if (!moveMesh_)
    {
        return;
    }

    const scalarField newV(mass0/rho_);

    Info<< "Initial/final volumes = " << gSum(regionMesh().V()) << ", "
        << gSum(newV) << " [m3]" << endl;

    // move the mesh
    const labelList moveMap = moveMesh(regionMesh().V() - newV, minimumDelta_);

    // flag any cells that have not moved as non-reacting
    forAll(moveMap, i)
    {
        if (moveMap[i] == 0)
        {
            solidChemistry_->setCellReacting(i, false);
        }
    }
}


void myReactingOneDim::solveContinuity()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    const scalarField mass0 = rho_*regionMesh().V();

    fvScalarMatrix rhoEqn
    (
        fvm::ddt(rho_) == -solidChemistry_->RRg()
    );

    if (regionMesh().moving())
    {
        surfaceScalarField phiRhoMesh
        (
            fvc::interpolate(rho_)*regionMesh().phi()
        );

        rhoEqn += fvc::div(phiRhoMesh);
    }

    rhoEqn.solve();

    updateMesh(mass0);
}


void myReactingOneDim::solveSpeciesMass()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    volScalarField Yt(0.0*Ys_[0]);

    for (label i=0; i<Ys_.size()-1; i++)
    {
        volScalarField& Yi = Ys_[i];

        fvScalarMatrix YiEqn
        (
            fvm::ddt(rho_, Yi) == solidChemistry_->RRs(i)
        );

        if (regionMesh().moving())
        {
            surfaceScalarField phiYiRhoMesh
            (
                fvc::interpolate(Yi*rho_)*regionMesh().phi()
            );

            YiEqn += fvc::div(phiYiRhoMesh);

        }

        YiEqn.solve("Yi");
        Yi.max(0.0);
        Yt += Yi;
    }

    Ys_[Ys_.size() - 1] = 1.0 - Yt;

}


void myReactingOneDim::solveEnergy()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    fv::options& fvOptions(fv::options::New(regionMesh()));

    dimensionedScalar Cp0("Cp0", dimEnergy/dimMass/dimTemperature, solidThermo_->composition().Cp(0, 1.0e05, 300.0));
    dimensionedScalar Cp1("Cp1", dimEnergy/dimMass/dimTemperature, solidThermo_->composition().Cp(1, 1.0e05, 300.0));   

    tmp<volScalarField> alpha(solidThermo_->alpha());

    fvScalarMatrix hEqn
    (
        fvm::ddt(rho_, h_)
      - fvm::laplacian(alpha, h_)
      + fvc::laplacian(alpha, h_)
      - fvc::laplacian(kappa(), T())
     ==
        chemistryQdot_
      + solidChemistry_->RRs(0)*T()*Cp0
      + solidChemistry_->RRs(1)*T()*Cp1
      + fvOptions(rho_, h_)
    );

    if (gasHSource_)
    {
        const surfaceScalarField phiGas(fvc::interpolate(phiHsGas_));
        hEqn += fvc::div(phiGas);
    }

    if (qrHSource_)
    {
        const surfaceScalarField phiqr(fvc::interpolate(qr_)*nMagSf());
        hEqn += fvc::div(phiqr);
    }

    if (charOxiSource_)
    {
        hEqn += charOxiQdot_;
    }

    if (regionMesh().moving())
    {
        surfaceScalarField phihMesh
        (
            fvc::interpolate(rho_*h_)*regionMesh().phi()
        );

        hEqn += fvc::div(phihMesh);
    }

    hEqn.relax();
    hEqn.solve();
}


void myReactingOneDim::calculateMassTransfer()
{
    totalGasMassFlux_ = 0;
    forAll(intCoupledPatchIDs_, i)
    {
        const label patchi = intCoupledPatchIDs_[i];
        totalGasMassFlux_ += gSum(phiGas_.boundaryField()[patchi]);
    }

    if (infoOutput_)
    {
        totalHeatRR_ = fvc::domainIntegrate(chemistryQdot_);

        addedGasMass_ +=
            fvc::domainIntegrate(solidChemistry_->RRg())*time_.deltaT();
        lostSolidMass_ +=
            fvc::domainIntegrate(solidChemistry_->RRs())*time_.deltaT();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

myReactingOneDim::myReactingOneDim
(
    const word& modelType,
    const fvMesh& mesh,
    const word& regionType
)
:
    pyrolysisModel(modelType, mesh, regionType),
    solidThermo_(solidReactionThermo::New(regionMesh())),
    solidChemistry_(basicSolidChemistryModel::New(solidThermo_())),
    radiation_(radiationModel::New(solidThermo_->T())),
    rho_
    (
        IOobject
        (
            "rho",
            regionMesh().time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        solidThermo_->rho()
    ),
    Ys_(solidThermo_->composition().Y()),
    h_(solidThermo_->he()),
    nNonOrthCorr_(-1),
    maxDiff_(10),
    minimumDelta_(1e-4),

    phiGas_
    (
        IOobject
        (
            "phiGas",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime, 0.0)
    ),

    phiHsGas_
    (
        IOobject
        (
            "phiHsGas",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime, 0.0)
    ),

    chemistryQdot_
    (
        IOobject
        (
            "chemistryQdot",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),

    qr_
    (
        IOobject
        (
            "qr",
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),

    Xchar_
    (
        IOobject
        (
            "Xchar",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),

    mChar_
    (
        IOobject
        (
            "mChar",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass, 0.0)
    ),

    mCharBurnt_
    (
        IOobject
        (
            "mCharBurnt",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass, 0.0)
    ),
    charOxiQdot_
    (
        IOobject
        (
            "charOxiQdot",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),

    phiO2_
    (
        IOobject
        (
            "phiO2",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime, 0.0)
    ),

    phiCO2_
    (
        IOobject
        (
            "phiCO2",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime, 0.0)
    ),

    lostSolidMass_(dimensionedScalar("zero", dimMass, 0.0)),
    addedGasMass_(dimensionedScalar("zero", dimMass, 0.0)),
    totalGasMassFlux_(0.0),
    totalHeatRR_(dimensionedScalar("zero", dimEnergy/dimTime, 0.0)),
    gasHSource_(false),
    qrHSource_(false),
    charOxiSource_(false),
    useChemistrySolvers_(true)
{
    if (active_)
    {
        read();
    }
}


myReactingOneDim::myReactingOneDim
(
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& regionType
)
:
    pyrolysisModel(modelType, mesh, dict, regionType),
    solidThermo_(solidReactionThermo::New(regionMesh())),
    solidChemistry_(basicSolidChemistryModel::New(solidThermo_())),
    radiation_(radiationModel::New(solidThermo_->T())),
    rho_
    (
        IOobject
        (
            "rho",
            regionMesh().time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        solidThermo_->rho()
    ),
    Ys_(solidThermo_->composition().Y()),
    h_(solidThermo_->he()),
    nNonOrthCorr_(-1),
    maxDiff_(10),
    minimumDelta_(1e-4),

    phiGas_
    (
        IOobject
        (
            "phiGas",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime, 0.0)
    ),

    phiHsGas_
    (
        IOobject
        (
            "phiHsGas",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime, 0.0)
    ),

    chemistryQdot_
    (
        IOobject
        (
            "chemistryQdot",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),

    qr_
    (
        IOobject
        (
            "qr",
            time().timeName(),
            regionMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),

    Xchar_
    (
        IOobject
        (
            "Xchar",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),

    mChar_
    (
        IOobject
        (
            "mChar",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass, 0.0)
    ),

    mCharBurnt_
    (
        IOobject
        (
            "mCharBurnt",
            time().timeName(),
            regionMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass, 0.0)
    ),
    charOxiQdot_
    (
        IOobject
        (
            "charOxiQdot",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimEnergy/dimTime/dimVolume, 0.0)
    ),

    phiO2_
    (
        IOobject
        (
            "phiO2",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime, 0.0)
    ),

    phiCO2_
    (
        IOobject
        (
            "phiCO2",
            time().timeName(),
            regionMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar("zero", dimMass/dimTime, 0.0)
    ),

    lostSolidMass_(dimensionedScalar("zero", dimMass, 0.0)),
    addedGasMass_(dimensionedScalar("zero", dimMass, 0.0)),
    totalGasMassFlux_(0.0),
    totalHeatRR_(dimensionedScalar("zero", dimEnergy/dimTime, 0.0)),
    gasHSource_(false),
    qrHSource_(false),
    charOxiSource_(false),
    useChemistrySolvers_(true)
{
    if (active_)
    {
        read(dict);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

myReactingOneDim::~myReactingOneDim()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar myReactingOneDim::addMassSources(const label patchi, const label facei)
{
    label index = 0;
    forAll(primaryPatchIDs_, i)
    {
        if (primaryPatchIDs_[i] == patchi)
        {
            index = i;
            break;
        }
    }

    const label localPatchId =  intCoupledPatchIDs_[index];

    const scalar massAdded = phiGas_.boundaryField()[localPatchId][facei];

    if (debug)
    {
        Info<< "\nPyrolysis region: " << type() << "added mass : "
            << massAdded << endl;
    }

    return massAdded;
}


scalar myReactingOneDim::solidRegionDiffNo() const
{
    scalar DiNum = -great;

    if (regionMesh().nInternalFaces() > 0)
    {
        surfaceScalarField KrhoCpbyDelta
        (
            sqr(regionMesh().surfaceInterpolation::deltaCoeffs())
           *fvc::interpolate(kappa())
           /fvc::interpolate(Cp()*rho_)
        );

        DiNum = max(KrhoCpbyDelta.primitiveField())*time().deltaTValue();
    }

    return DiNum;
}


scalar myReactingOneDim::maxDiff() const
{
    return maxDiff_;
}


const volScalarField& myReactingOneDim::rho() const
{
    return rho_;
}


const volScalarField& myReactingOneDim::T() const
{
    return solidThermo_->T();
}


const tmp<volScalarField> myReactingOneDim::Cp() const
{
    return solidThermo_->Cp();
}


tmp<volScalarField> myReactingOneDim::kappaRad() const
{
    return radiation_->absorptionEmission().a();
}


tmp<volScalarField> myReactingOneDim::kappa() const
{
    return solidThermo_->kappa();
}


const surfaceScalarField& myReactingOneDim::phiGas() const
{
    return phiGas_;
}


void myReactingOneDim::preEvolveRegion()
{
    pyrolysisModel::preEvolveRegion();

    // Initialise all cells as able to react
    forAll(h_, celli)
    {
        solidChemistry_->setCellReacting(celli, true);
    }
}


void myReactingOneDim::evolveRegion()
{
    Info<< "\nEvolving pyrolysis in region: " << regionMesh().name() << endl;

    if (useChemistrySolvers_)
    {
        solidChemistry_->solve(time().deltaTValue());
    }
    else
    {
        solidChemistry_->calculate();
    }

    solveContinuity();

    chemistryQdot_ = solidChemistry_->Qdot()();

    updateFields();

    solveSpeciesMass();

    for (int nonOrth=0; nonOrth<=nNonOrthCorr_; nonOrth++)
    {
        solveEnergy();
    }

    calculateMassTransfer();

    solidThermo_->correct();

    Info<< "pyrolysis min/max(T) = "
        << gMin(solidThermo_->T().primitiveField())
        << ", "
        << gMax(solidThermo_->T().primitiveField())
        << endl;
}


void myReactingOneDim::info()
{
    Info<< "\nPyrolysis in region: " << regionMesh().name() << endl;

    Info<< indent << "Total gas mass produced  [kg] = "
        << addedGasMass_.value() << nl
        << indent << "Total solid mass lost    [kg] = "
        << lostSolidMass_.value() << nl
        << indent << "Total pyrolysis gases  [kg/s] = "
        << totalGasMassFlux_ << nl
        << indent << "Total heat release rate [J/s] = "
        << totalHeatRR_.value() << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace regionModels
} // End namespace pyrolysisModels

// ************************************************************************* //
