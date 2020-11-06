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

\*---------------------------------------------------------------------------*/

#include "multiStep.H"
//#include "singleStepReactingMixture.H"
#include "reactingMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvcDiv.H"
#include "fvmLaplacian.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


template<class ThermoType>
const Foam::reactingMixture<ThermoType>&
Foam::radiationModels::sootModels::multiStep<ThermoType>::checkThermo
(
    const fluidThermo& thermo
)
{
    if (isA<reactingMixture<ThermoType>>(thermo))
    {
        return dynamic_cast<const reactingMixture<ThermoType>&>
        (
            thermo
        );
    }
    else
    {
        FatalErrorInFunction
            << "Inconsistent thermo package for " << thermo.type()
            << "Please select a thermo package based on "
            << "reactingMixture" << exit(FatalError);

        return dynamic_cast<const reactingMixture<ThermoType>&>
        (
            thermo
        );
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::radiationModels::sootModels::multiStep<ThermoType>::multiStep
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& modelType
)
:
    sootModel(dict, mesh, modelType),

    soot_
    (
        IOobject
        (
            "soot",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    PDSoot_
    (
        IOobject
        (
            "PDSoot",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    dNSootdt_
    (
        IOobject
        (
            "dNSootdt",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("dNSootdt", dimensionSet(0,-3,-1,0,0,0,0), 0.0)
    ),

    SootVF_
    (
        IOobject
        (
            "SootVF",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("dNSootdt", dimensionSet(0,-3,-1,0,0,0,0), 0.0)
    ),

    SSoot_
    (
        IOobject
        (
            "SSoot",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("SSoot", dimensionSet(0,-1,0,0,0,0,0), 0.0)
    ),

    dMSootdt_
    (
        IOobject
        (
            "dMSootdt",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("dMSootdt", dimensionSet(1,-3,-1,0,0,0,0), 0.0)
    ),

    r1_
    (
        IOobject
        (
            "r1",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("r1", dimensionSet(0,-3,-1,0,1,0,0), 0.0)
    ),

    r2_
    (
        IOobject
        (
            "r2",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("r2",dimensionSet(0,-3,-1,0,0,0,0), 0.0)
    ),

    r3_
    (
        IOobject
        (
            "r3",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("r3", dimensionSet(0,-3,-1,0,1,0,0), 0.0)
    ),

    r4_
    (
        IOobject
        (
            "r4",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("r4", dimensionSet(0,-3,-1,0,1,0,0), 0.0)
    ),

    r5_
    (
        IOobject
        (
            "r5",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("r5", dimensionSet(0,-3,-1,0,1,0,0), 0.0)
    ),

    MSoot_
    (
        IOobject
        (
            "MSoot",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("MSoot", dimensionSet(1,-3,0,0,0,0,0), 0.0)
    ),

    NSoot_
    (
        IOobject
        (
            "NSoot",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("NSoot", dimensionSet(0,-3,0,0,0,0,0), 0.0)
    ),

    YPrec_
    (
        IOobject
        (
            "YPrec",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("YPrec", dimensionSet(0,0,0,0,0,0,0), 0.0)
    ),

    coeffsDict_(dict.subOrEmptyDict(modelType + "Coeffs")),
    Sct_(readScalar(coeffsDict_.lookup("Sct"))),
    MWc_(readScalar(coeffsDict_.lookup("MWc"))),
    wPrec_(readScalar(coeffsDict_.lookup("wPrec"))),
    rhoS_(readScalar(coeffsDict_.lookup("rhoS"))),
    Ca_(readScalar(coeffsDict_.lookup("Ca"))),
    Ta_(readScalar(coeffsDict_.lookup("Ta"))),
    Cb_(readScalar(coeffsDict_.lookup("Cb"))),
    Cg_(readScalar(coeffsDict_.lookup("Cg"))),
    Tg_(readScalar(coeffsDict_.lookup("Tg"))),
    m_(readScalar(coeffsDict_.lookup("m"))),
    q_(readScalar(coeffsDict_.lookup("q"))),
    Cw1_(readScalar(coeffsDict_.lookup("Cw1"))),
    CollEff_(readScalar(coeffsDict_.lookup("CollEff"))),
    Cw2_(readScalar(coeffsDict_.lookup("Cw2"))),
    Tw2_(readScalar(coeffsDict_.lookup("Tw2"))),
    NA_(readScalar(coeffsDict_.lookup("NA"))),
    thermo_(mesh.lookupObject<fluidThermo>(basicThermo::dictName))

{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::radiationModels::sootModels::multiStep<ThermoType>::~multiStep()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::radiationModels::sootModels::multiStep<ThermoType>::correct()
{
    //Look for the phi and rho
    const surfaceScalarField& phi_(mesh_.lookupObject<surfaceScalarField>("phi"));
    const volScalarField& rho_(mesh_.lookupObject<volScalarField>("rho"));

    //Read the rho T p from thermo 
    const volScalarField& T = thermo_.T();
    const volScalarField& p = thermo_.p();

    //Species/soot properties, reference conditions and other constants
    scalar kb = 1.3806488e-23;//????? some article is e-16 
    scalar Pref = 101325;// reference pressure
    scalar wC2H2=26.03824;//relative molecular mass of C2H2
    scalar wOH=17.00737; //relative molecular mass of OH
    scalar wO2=31.9988;// this is relative molecular mass

    //Get index of C2H2 in mixture 
    const volScalarField& Y_C2H2 = mesh_.lookupObject<volScalarField>("C2H2");

    //Get index of O2 in mixture
    const volScalarField& Y_O2 = mesh_.lookupObject<volScalarField>("O2");

    forAll(T, celli)
    {
	scalar Ti = T[celli];
	
	if (soot_[celli]<1.0e-12)
        {		
            MSoot_[celli] = 0.0;
        }
        else
        {
            MSoot_[celli] = soot_[celli]*rho_[celli]; //soot mass density [kg/m3]
        }
	
        if (PDSoot_[celli]<1.0e-15)
        {		
            NSoot_[celli] = 0.0;
        }
	else
        {
            NSoot_[celli] = PDSoot_[celli]*rho_[celli]*NA_; //soot number density [1/m3]  
        }
		
	SootVF_[celli] = soot_[celli]*rho_[celli]/rhoS_; //soot volume fraction [-]
	SSoot_[celli] = pow(6*MSoot_[celli]/rhoS_,2.0/3.0)*pow(constant::mathematical::pi*NSoot_[celli],1.0/3.0); //specific surface area [1/m];

	r1_[celli] = Ca_*(Y_C2H2[celli]*rho_[celli]/wPrec_)*exp(-Ta_/Ti); //nucleation [kmol/m3/s]
	r2_[celli] = Cb_*sqrt(24*kb*T[celli]/rhoS_)*pow(6*MSoot_[celli]/constant::mathematical::pi/rhoS_,1.0/6.0)*pow(NSoot_[celli],11.0/6.0); //coagulation [1/m3/s]  

	if (m_ == 0.5)
        {		
            r3_[celli] = pow(p[celli]/Pref,q_)*Cg_*(Y_C2H2[celli]*rho_[celli]/wC2H2)*exp(-Tg_/Ti)*sqrt(SSoot_[celli]); //Surface growth
        }
        else if (m_ == 1.0)
        {
            r3_[celli] = pow(p[celli]/Pref,q_)*Cg_*(Y_C2H2[celli]*rho_[celli]/wC2H2)*exp(-Tg_/Ti)*SSoot_[celli]; //Surface growth
        }
        else
        {
            FatalErrorIn
            (
                ""
            )
            << "Invalid integer value for m" << nl
            << "Valid integer values are 0.5 and 1.0" << nl
            << "m = 0.5 indicates surface growth rate is a square root function of SSoot" << nl
            << "m = 1.0 indicates surface growth rate is a linear function SSoot" << nl
            << exit(FatalError);
	}

        if (Cw1_ > 0)
        {
            const volScalarField& Y_OH = mesh_.lookupObject<volScalarField>("OH");
	    r4_[celli] = Cw1_*CollEff_*(Y_OH[celli]*rho_[celli]/wOH)*sqrt(Ti)*SSoot_[celli]; //Oxidation due to OH [kmol/m3/s] 
        }
	
	r5_[celli] = Cw2_*(Y_O2[celli]*rho_[celli]/wO2)*exp(-Tw2_/Ti)*sqrt(Ti)*SSoot_[celli]; //Oxidation due to O2 [kmol/m3/s]	

	//Calculation for rates	
	dNSootdt_[celli] = NA_*r1_[celli]-r2_[celli];
	dMSootdt_[celli] = MWc_*(100*r1_[celli]+2*r3_[celli]-r4_[celli]-r5_[celli]);

	// Stabilizing model - the prevention of negative concentration
        scalar deltaT =0.0;
        deltaT = this->mesh().time().deltaTValue();
	
        if (MSoot_[celli]+dMSootdt_[celli]*deltaT<0.0)
            dMSootdt_[celli]=-MSoot_[celli]/deltaT;

        if (NSoot_[celli]+dNSootdt_[celli]*deltaT<0.0)
            dNSootdt_[celli]=-NSoot_[celli]/deltaT;
    }

    const compressible::turbulenceModel& 
    turbModel = rho_.db().lookupObject<compressible::turbulenceModel>
    (
        turbulenceModel::propertiesName
    );

    // Transport equation for particle number
    tmp<fvScalarMatrix> PDSootEqn 
    (
        fvm::ddt(rho_, PDSoot_)
      + fvm::div(phi_, PDSoot_)
      - fvm::laplacian(turbModel.muEff()/Sct_, PDSoot_)
      ==
        dNSootdt_/(NA_*dimensionedScalar(dimless/dimMoles, 1.0))
    );

    solve(PDSootEqn);

    // Transport equation for soot mass density
    tmp<fvScalarMatrix> sootEqn
    (
        fvm::ddt(rho_, soot_)
      + fvm::div(phi_, soot_)
      - fvm::laplacian(turbModel.muEff()/Sct_, soot_)
      ==
        dMSootdt_
    );

    solve(sootEqn);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
