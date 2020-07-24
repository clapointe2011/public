#include "costFunctionLibrary.H"
#include "fvCFD.H"

CostFunction::CostFunction(Foam::fvMesh &mesh)
    :
      mesh(mesh)
{}

Foam::scalar CostFunction::eval()
{
    const Foam::wordList costFunctionPatches(
            mesh.solutionDict().subDict("OPTIMIZATION").lookup("costFunctionPatches")
        );

    const Foam::string costFunction =
            mesh.solutionDict().subDict("OPTIMIZATION").lookup("costFunction");

    Foam::scalar J = 0;

    //Foam::Info << "Cost Function: " << costFunction << nl << Foam::endl;

    if(costFunction=="sonicFlowPoint"){
        J = eval_sonicFlowPoint(costFunctionPatches);
    }
    else if(costFunction=="incompressiblePressureDrop"){
        J = eval_incompPressureDrop(costFunctionPatches);
    }
    else if(costFunction=="sonicTotalPressure"){
	    J = eval_sonicTotalPressure(costFunctionPatches);
    }
    else if(costFunction=="rhoCentral"){
	    J = eval_rhoCentral(costFunctionPatches);
    }
    else if(costFunction=="reacSonicTotalPressure"){
	    J = eval_reacSonicTotalPressure(costFunctionPatches);
    }
    else if(costFunction=="minMaxFlowOutlet"){
	    J = eval_flowOutlet(costFunctionPatches);
    }
    else if(costFunction=="reacFlow"){
	    J = eval_reacFlow(costFunctionPatches);
    }
    else if(costFunction=="minMaxFlowSection"){
	    J = eval_flowSection(costFunctionPatches);
    }
    else if(costFunction=="minRollerTemp"){
	    J = eval_rollerTemp(costFunctionPatches);
    }
    else if(costFunction=="shockAngle"){
	    J = eval_shockAngle(costFunctionPatches);
    }
    else if(costFunction=="vorticity"){
	    J = eval_vorticity(costFunctionPatches);
    }
    else{
        Foam::Info << "Unknown Cost Function!" << Foam::endl;
    }

    return J;
}

Foam::scalar CostFunction::eval_incompPressureDrop(const Foam::wordList &costFunctionPatches)
{
    //NOTE: this is only pressure drop if evaluated for inlet and outlet of domain

    const Foam::volScalarField& p = mesh.lookupObject<Foam::volScalarField>("p");
    const Foam::volScalarField& alpha = mesh.lookupObject<Foam::volScalarField>("alpha");
    const Foam::surfaceScalarField& phi = mesh.lookupObject<Foam::surfaceScalarField>("phi");

    scalar J = 0;

    scalar F1 = 0;
    scalar F2 = 0;

    forAll(costFunctionPatches,cI)
    {
        Foam::label patchI = mesh.boundaryMesh().findPatchID(costFunctionPatches[cI]);
        const Foam::fvPatch& patch = mesh.boundary()[patchI];

	if (cI == 0){
	    F1 += gSum
		(
		    - phi.boundaryField()[patchI]*(p.boundaryField()[patchI] + 0.5*magSqr(phi.boundaryField()[patchI]/patch.magSf()))
                );
	}else{
	    F2 += gSum
		(
		    - phi.boundaryField()[patchI]*(p.boundaryField()[patchI] + 0.5*magSqr(phi.boundaryField()[patchI]/patch.magSf()))
                );
	}
    }

    scalar alphaMax = 1e2;
    scalar alphaTot = gSum(alpha);
    scalar numCell = mesh.cells().size();
    scalar maxAlphaTot = alphaMax*numCell;
    scalar Vtar = 0.04;

    J = pow((1e2*F1-1e2*F2),2) + 0*pow(((1-(alphaTot/maxAlphaTot))-Vtar),2);

    return J;
}

Foam::scalar CostFunction::eval_shockAngle(const Foam::wordList &costFunctionPatches)
{
    const Foam::volVectorField& U = mesh.lookupObject<Foam::volVectorField>("U");
    const Foam::volScalarField& alpha = mesh.lookupObject<Foam::volScalarField>("alpha");
    const Foam::volScalarField& p = mesh.lookupObject<Foam::volScalarField>("p");
    const Foam::surfaceScalarField& Maf = mesh.lookupObject<Foam::surfaceScalarField>("((mag(phi)|rhof)|((cSf_pos*a_pos)+(cSf_neg*a_neg)))");

    const volScalarField& Ma = fvc::average(Maf);

    //const volVectorField& gradP = Foam::fvc::grad(p);

    //initialize variables

    scalar J = 0;
    //scalar J1 = 0;
    //scalar J2 = 0;

    //scalar UxWave = 0;
    //scalar UyWave = 0;
    //scalar UxP = 0;
    //scalar UyP = 0;
    scalar Uin = 560;

    scalar pIn = 23000;
    scalar wavePRatio = 0;

    //scalar gamma = 1.4;
    scalar waveMa = 0;
    scalar waveAngle = 0;
    //scalar waveUfluc = 0;
    //scalar alphaPgrad = 0;
    //scalar alphaP = 0;
    //scalar alphaX = 0;
    //scalar alphaY = 0;

    scalar xPos = 0;
    scalar yPos = 0;
    scalar alphaPos = 0;
    scalar UxPos = 0;
    scalar UyPos = 0;
    //scalar UmagPos = 0;
    scalar pPos = 0;
    //scalar rhoPos = 0;
    //scalar aPos = 0;
    scalar MaPos = 0;
    //scalar Umean = 0;
    scalar n = SMALL;
    //scalar m = SMALL;

    dimensionedScalar alphaMax(mesh.solutionDict().subDict("OPTIMIZATION").lookup("alphaMax")); //CL
    
    //read targets 

    scalar targetMa = readScalar(mesh.solutionDict().subDict("OPTIMIZATION").lookup("targetMa")); //CL
    //scalar targetVolFrac = readScalar(mesh.solutionDict().subDict("OPTIMIZATION").lookup("targetVolFrac")); //CL
    scalar targetAngle = readScalar(mesh.solutionDict().subDict("OPTIMIZATION").lookup("targetAngle")); //CL
    //scalar targetPRatio = readScalar(mesh.solutionDict().subDict("OPTIMIZATION").lookup("targetPRatio")); //CL

    const Foam::vector positionStart = Foam::vector(
	    mesh.solutionDict().subDict("OPTIMIZATION").lookup("rampStart")
	);

    scalar xPosRampStart = positionStart.component(0);
    scalar yPosRampStart = positionStart.component(1);

    /*const Foam::vector direction = Foam::vector(
            mesh.solutionDict().subDict("OPTIMIZATION").lookup("flowDirection")
        );*/

    //scalar maxYposAlpha = 0;
    //scalar maxXposAlpha = 0;
    //label maxCellxI;
    //label maxCellyI;
    //scalar maxAlpha = 0;
    //label maxCellI;

    /*forAll(alpha, cellI){
	if (alpha[cellI] > 0){
	    if (maxYposAlpha < mesh.C()[cellI].component(1)){
		maxYposAlpha = mesh.C()[cellI].component(1);
	    }
	}
    }*/

    //scalar maxAlphaXpos = mesh.C()[maxCellI].component(0);
    //scalar maxAlphaYpos = mesh.C()[maxCellI].component(1);

    //scalar thetaMaxAlpha = Foam::atan(maxAlphaYpos/maxAlphaXpos)*180/Foam::constant::mathematical::pi;

    /*forAll(U,i){
	xPos = mesh.C()[i].component(0);
	yPos = mesh.C()[i].component(1);
	alphaPos = alpha.internalField()[i];
	UxPos = U.internalField()[i].component(0);

	if (xPos > xPosRampStart && yPos > yPosRampStart){
	    if (alphaPos == 0){
		if (UxPos < Uin){
		    Umean += mag(U.internalField()[i]);
		    m = m + 1;
		}
	    }
	}
    }

    Umean = Umean/m;*/

    forAll(costFunctionPatches,pI){

        /*Foam::label patchI = mesh.boundaryMesh().findPatchID(costFunctionPatches[pI]);
        const Foam::fvPatch& patch = mesh.boundary()[patchI];

	    forAll(U.boundaryField()[patchI],pI){
	        if(alpha.boundaryField()[patchI][pI] == 0){
	        	UxP += U.boundaryField()[patchI][pI].component(0);
	        	UyP += U.boundaryField()[patchI][pI].component(1);
		    nP = nP + 1;
	        }
	    }*/

	    forAll(U,i){

	        /*const labelList& neighbourCells = mesh.cellCells()[i];
	        bool enclosed = true;

	        forAll(neighbourCells,nI){
		    if (alpha[nI] == 0){
		        enclosed = false;
		        break;
		    }
	        }*/

	        xPos = mesh.C()[i].component(0);
	        yPos = mesh.C()[i].component(1);
	        alphaPos = alpha.internalField()[i];
	        UxPos = U.internalField()[i].component(0);
	        UyPos = U.internalField()[i].component(1);
	        pPos = p.internalField()[i];
	        MaPos = Ma[i];

	        /*if (xPos > xPosRampStart && yPos > yPosRampStart){
		    if (alphaPos == 0 && !enclosed){
		        if (UxPos < Uin){
			    waveAngle += Foam::atan(UyPos/UxPos)*180/Foam::constant::mathematical::pi;
			    waveMa += MaPos;
			    waveUfluc += UmagPos;
			    n = n + 1;
		        }
		    }
	        }*/

	        if (xPos > xPosRampStart && yPos > yPosRampStart){
		        if (alphaPos == 0){
		            if (UxPos < Uin){
			            waveAngle += Foam::atan(UyPos/UxPos)*180/Foam::constant::mathematical::pi;
			            waveMa += MaPos;
			            wavePRatio += pPos/pIn;
			            n = n + 1;
		            }
		        }
	        }
	    }
    }

    //calculate J

    /*// read weights
    scalar weightMa = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightMa",1); //CL
    scalar weightAngle = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightAngle",1); //CL
    scalar weightUfluc = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightUfluc",1); //CL
    scalar weightVolFrac = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightVolFrac",1); //CL
    //scalar weightPgrad = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightPgrad",1); //CL
    //scalar weightP = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightP",1); //CL
    //scalar weightPacking = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightPacking",1); //CL

    J =   weightMa*pow(waveMa-targetMa*n,2)
	+ weightAngle*pow(waveAngle-targetAngle*n,2)
	+ weightUfluc*pow(waveUfluc-Umean*n,2)
	//+ weightPgrad*alphaPgrad
	//+ weightP*alphaP
	//+ weightPacking*alphaY
	+ weightVolFrac*pow((gSum(alpha)/(alphaMax.value()*mesh.cells().size()))-targetVolFrac,2);*/

    // read weights
    scalar relativeWeight = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("relativeWeight",0.5);

    J =   relativeWeight*pow((waveMa/n - targetMa),2)
	+ (1-relativeWeight)*pow((waveAngle/n - targetAngle),2);

    /*J =   relativeWeight*pow((wavePRatio/n - targetPRatio),2)
	+ (1-relativeWeight)*pow((waveAngle/n - targetAngle),2);*/


    Info << "\nmean waveAngle : " << waveAngle/n << " mean waveMa : " << waveMa/n <<  " mean PRatio : " << wavePRatio/n << nl << endl;

    return J;
}

Foam::scalar CostFunction::eval_rollerTemp(const Foam::wordList &costFunctionPatches)
{
    const Foam::volScalarField& T = mesh.lookupObject<Foam::volScalarField>("T");

    scalar J = 0;

    scalar intTargetTemp = 0;
    scalar intTemp = 0;

    forAll(costFunctionPatches,pI)
    {

        Foam::label patchI = mesh.boundaryMesh().findPatchID(costFunctionPatches[pI]);
        const Foam::fvPatch& currPatch = mesh.boundary()[patchI];

	/*forAll(currPatch,facei){
	    Foam::label faceCelli = currPatch.faceCells()[facei];
	    J += pow(T[facei]-targetTemp,2);
	}*/

	intTemp += gSum(T.boundaryField()[patchI]);
	intTargetTemp += 500*currPatch.size();
    }

    J = pow(intTemp - intTargetTemp,2);

    return J;
}

Foam::scalar CostFunction::eval_flowOutlet(const Foam::wordList &costFunctionPatches)
{
    const Foam::volVectorField& U = mesh.lookupObject<Foam::volVectorField>("U");
    const Foam::volScalarField& alpha = mesh.lookupObject<Foam::volScalarField>("alpha");
    const Foam::volScalarField& p = mesh.lookupObject<Foam::volScalarField>("p");

    scalar J = 0;

    volTensorField gradU = Foam::fvc::grad(U);
    //volScalarField gradUx = gradU.component(0);
    //volScalarField gradUy = gradU.component(1);

    //volVectorField gradp = Foam::fvc::grad(p);
    //volScalarField gradpx = gradp.component(0);
    //volScalarField gradpy = gradp.component(1);

    //volVectorField dAdx = Foam::fvc::grad(alpha);
    //volScalarField dAdx_x = dAdx.component(0);
    //volScalarField dAdx_y = dAdx.component(1);

    //volTensorField dA2dx2 = Foam::fvc::grad(dAdx);
    //volScalarField dA2dx2_x = dA2dx2.component(0);
    //volScalarField dA2dx2_y = dA2dx2.component(4);


    scalar cutOffOutletV = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("cutOffOutletV",0);
    //scalar targetAngle = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("targetAngle",0);

    scalar recircU = 0;
    scalar outletLV = 0;
    //scalar alphaCellCount = 0;
    //scalar internalV = 0;
    //scalar sumOutletAngle = 0;
    //scalar sumTargetAngle = 0;
    //scalar meanOutletAngle = 0;
    scalar gradientU = 0;
    //scalar gradientP = 0;
    //scalar gradientA = 0;
    //scalar curvatureA = 0;

    forAll(costFunctionPatches,cI)
    {
        Foam::label patchI = mesh.boundaryMesh().findPatchID(costFunctionPatches[cI]);
        const Foam::fvPatch& currPatch = mesh.boundary()[patchI];

	forAll(currPatch,facei){
	    Foam::label faceCelli = currPatch.faceCells()[facei];
	    Foam::vector tmp = U[faceCelli];

	    /*if (alpha[faceCelli] == 0){
	        sumOutletAngle += Foam::atan(tmp.component(1)/tmp.component(0))*180/Foam::constant::mathematical::pi;
	    }*/

	    if (tmp.component(1) < cutOffOutletV){
		outletLV += mag(tmp.component(1));
	    }
	}
	
	//sumTargetAngle = targetAngle*currPatch.size();
	//meanOutletAngle = sumOutletAngle/currPatch.size();
    }

    forAll(mesh.C(),celli){

	if (U[celli].component(0) < 0){
	    recircU += mag(U[celli].component(0));
	}//else if (alpha[celli] != 0){
	//    alphaCellCount += 1;
	//}

	gradientU += (mag(gradU[celli].component(0)) + mag(gradU[celli].component(1)));
	//gradientP += (mag(gradpx[celli]) + mag(gradpy[celli]));
	//gradientA += (mag(gradAlphax[celli]) + mag(gradAlphay[celli]));
	
	//if(alpha[celli] != 0){
	//curvatureA += (mag(dA2dx2_x[celli]) + mag(dA2dx2_y[celli]));
	//}
    }

    //scalar targetAlphaFrac = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("targetAlphaFrac",0.3);
    //dimensionedScalar alphaMax(mesh.solutionDict().subDict("OPTIMIZATION").lookup("alphaMax")); //CL
    //scalar alphaFrac = gSum(alpha)/(alphaMax.value()*alpha.size());
    //scalar alphaFrac = alphaCellCount/alpha.size();

    //Info << "\nmagnitude recirculation : " << recircU << nl << endl;
    //Info << "magnitude above cutoff outlet : " << outletHU << nl << endl;
    //Info << "magnitude V below cutoff outlet : " << outletLV << nl << endl;
    //Info << "magnitude V below cutoff  internal : " << internalV << nl << endl;
    //Info << "mean angle at outlet: " << meanOutletAngle << nl << endl;
    //Info << "volume fraction filled : " << alphaTot/maxAlphaTot << nl << endl;

    // read weights
    scalar weightRecirc = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightRecirc",1); //CL
    scalar weightGradU = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightGradU",0); //CL
    //scalar weightAlphaMag = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightAlphaMag",0); //CL
    //scalar weightAlphaFrac = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightAlphaFrac",0); //CL
    //scalar weightGradP = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightGradP",0); //CL
    scalar weightOutletLV = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightOutletLV",0); //CL
    //scalar weightOutletAngle = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightOutletAngle",0); //CL
    //scalar weightGradAlpha = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightGradAlpha",0); //CL
    //scalar weightCurvatureA = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightCurvatureA",0); //CL

    J =   weightRecirc*recircU
	+ weightOutletLV*outletLV
	//+ weightOutletAngle*abs(sumOutletAngle-sumTargetAngle)
	+ weightGradU*gradientU;
	//+ weightAlphaMag*gSum(alpha)
	//+ weightAlphaFrac*pow((alphaFrac - targetAlphaFrac),2)
	//+ weightGradP*gradientP;
	//+ weightGradAlpha*gradientA
	//+ weightCurvatureA*curvatureA;

    return J;
}

Foam::scalar CostFunction::eval_reacFlow(const Foam::wordList &costFunctionPatches)
{
    const Foam::volVectorField& U = mesh.lookupObject<Foam::volVectorField>("U");
    const Foam::volScalarField& H2O = mesh.lookupObject<Foam::volScalarField>("H2O");
    const Foam::volScalarField& CO2 = mesh.lookupObject<Foam::volScalarField>("CO2");
    const Foam::volScalarField& CH4 = mesh.lookupObject<Foam::volScalarField>("CH4");
    const Foam::volScalarField& T = mesh.lookupObject<Foam::volScalarField>("T");
    const Foam::volScalarField& psi = mesh.lookupObject<Foam::volScalarField>("thermo:psi");
    const Foam::volScalarField& gamma = mesh.lookupObject<Foam::volScalarField>("gamma");
    const Foam::volScalarField& p = mesh.lookupObject<Foam::volScalarField>("p");
    const Foam::surfaceScalarField& phi = mesh.lookupObject<Foam::surfaceScalarField>("phi");

    scalar weightTotalPressure = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightTotalPressure",1); //CL
    scalar weightTemperature = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightTemperature",1); //CL

    scalar J = 0;

    //scalar recircU = 0;
    //scalar outletProd = 0;

    forAll(costFunctionPatches,cI)
    {
        Foam::label patchI = mesh.boundaryMesh().findPatchID(costFunctionPatches[cI]);
        const Foam::fvPatch& currPatch = mesh.boundary()[patchI];

	/*forAll(currPatch,facei){
	    Foam::label faceCelli = currPatch.faceCells()[facei];

	    if (alpha[faceCelli] == 0){
	        sumOutletAngle += Foam::atan(tmp.component(1)/tmp.component(0))*180/Foam::constant::mathematical::pi;
	    }

	    J -= gSum(
		p.boundaryField()[patchI]*pow(
		    (1+(gamma[faceCelli]-1/(2*gamma[faceCelli]))*psi[faceCelli]*magSqr(U[faceCelli])),
		    (gamma[faceCelli]/(gamma[faceCelli]-1))
		)
	    );

	}*/

	J -= weightTotalPressure*gSum(
		p.boundaryField()[patchI]*pow(
		    (1+(gamma.boundaryField()[patchI]-1/(2*gamma.boundaryField()[patchI]))*psi.boundaryField()[patchI]*magSqr(U.boundaryField()[patchI])),
		    (gamma.boundaryField()[patchI]/(gamma.boundaryField()[patchI]-1))
		)
	    );

	//J -= gSum(H2O.boundaryField()[patchI] + CO2.boundaryField()[patchI]);

	J += weightTemperature*gSum(phi.boundaryField()[patchI]*T.boundaryField()[patchI]);
	
    }

    //J += gSum(CH4.internalField());

    /*forAll(mesh.C(),celli){

	if (U[celli].component(0) < 0){
	    J += mag(U[celli].component(0));
	}
    }*/

    /*forAll(mesh.C(),celli){
	Foam::vector tmp = U.internalField()[celli];

	if (tmp.component(0) < 0){
	    recircU += mag(tmp);
	}
    }*/

    //Info << "\nmagnitude recirculation : " << recircU << nl << endl;

    // read weights
    //scalar weightRecirc = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightRecirc",1); //CL
    //scalar weightProd = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightProd",1); //CL

    //J = weightRecirc*recircU
    //  - weightProd*outletProd;

    //J = weightRecirc*recircU;

    return J;
}

Foam::scalar CostFunction::eval_flowSection(const Foam::wordList &costFunctionPatches)
{
    const Foam::volVectorField& U = mesh.lookupObject<Foam::volVectorField>("U");

    scalar J = 0;

    const Foam::vector p1 = Foam::vector(
            mesh.solutionDict().subDict("OPTIMIZATION").lookup("bottomLeft")
        );

    const Foam::vector p2 = Foam::vector(
            mesh.solutionDict().subDict("OPTIMIZATION").lookup("topRight")
        );

    volTensorField gradU = Foam::fvc::grad(U);
    volScalarField gradUx = gradU.component(0);
    volScalarField gradUy = gradU.component(1);

    // read target velocity
    scalar cutOffU = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("cutOffU",0);

    scalar recircU = 0;
    scalar sectionL = 0;
    scalar sectionH = 0;
    scalar flowGradU = 0;
    scalar n = 0;
    //scalar outletAngle = 0;

    forAll(U,cI){
	scalar xPos = mesh.C()[cI].component(0);
	scalar yPos = mesh.C()[cI].component(1);

	if (U.internalField()[cI].component(0) < 0){
	    recircU += mag(U[cI]);
	}

	if (p1.component(0) <= xPos && xPos <= p2.component(0)){
	    if (p1.component(1) <= yPos && yPos <= p2.component(1)){
		if (U.internalField()[cI].component(1) < cutOffU){
		    sectionL += mag(U[cI]);
		    n = n + 1;
		}
	    }
	}
    }

    //Info << "\nrecirculation magnitude : " << recircU << " mean angle at outlet : " << outletAngle << nl << endl;

    // read weights
    scalar weightRecirc = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightRecirc",1); //CL
    scalar weightSectionL = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightSectionL",1); //CL
    //scalar weightAngle = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightAngle",1); //CL

    J = weightRecirc*recircU
      - weightSectionL*sectionL;

    return J;
}

Foam::scalar CostFunction::eval_rhoCentral(const Foam::wordList &costFunctionPatches)
{
    const Foam::volScalarField& p = mesh.lookupObject<Foam::volScalarField>("p");
    const Foam::surfaceScalarField& phi = mesh.lookupObject<Foam::surfaceScalarField>("phi");
    const Foam::volScalarField& psi = mesh.lookupObject<Foam::volScalarField>("thermo:psi");
    const Foam::volVectorField& U = mesh.lookupObject<Foam::volVectorField>("U");

    scalar J = 0;
    scalar pIn = 100000;
    scalar pOut = 0;
    scalar gamma = 1.4;

    /*forAll(costFunctionPatches,cI)
    {
        Foam::label patchI = mesh.boundaryMesh().findPatchID(costFunctionPatches[cI]);
        const Foam::fvPatch& currPatch = mesh.boundary()[patchI];

	    if (costFunctionPatches[cI] == "inlet"){
	        pIn += gSum(p.boundaryField()[patchI] + 0.5*magSqr(U.boundaryField()[patchI]));
	    }else{
	        pOut += gSum(p.boundaryField()[patchI]
	                    * pow((1+((gamma-1)/(2*gamma))*psi.boundaryField()[patchI]*magSqr(U.boundaryField()[patchI])),(gamma/(gamma-1))));
	    }
    }
    
    J += pow(pIn - pOut, 2);*/
    
    forAll(costFunctionPatches,cI)
    {
        Foam::label patchI = mesh.boundaryMesh().findPatchID(costFunctionPatches[cI]);
        const Foam::fvPatch& currPatch = mesh.boundary()[patchI];

        J += pow(pIn - gSum(p.boundaryField()[patchI]),2);
    }

    return J;
}

Foam::scalar CostFunction::eval_sonicTotalPressure(const Foam::wordList &costFunctionPatches)
{
    const Foam::volScalarField& p = mesh.lookupObject<Foam::volScalarField>("p");
    const Foam::surfaceScalarField& phi = mesh.lookupObject<Foam::surfaceScalarField>("phi");
    const Foam::volScalarField& psi = mesh.lookupObject<Foam::volScalarField>("thermo:psi");
    const Foam::volVectorField& U = mesh.lookupObject<Foam::volVectorField>("U");
    //const Foam::volScalarField& rho = mesh.lookupObject<Foam::volScalarField>("rho");
    //const Foam::surfaceScalarField& Maf = mesh.lookupObject<Foam::surfaceScalarField>("((mag(phi)|rhof)|((cSf_pos*a_pos)+(cSf_neg*a_neg)))");
    //const Foam::volScalarField& cSound = mesh.lookupObject<Foam::volScalarField>("cSound");

    scalar J = 0;
    //scalar outletTotalPressure = 0;
    scalar gamma = 1.4;
    //scalar outletP = 0;
    //scalar inletP = 0;
    //scalar outletLV = 0;

    //scalar cutOffOutletV = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("cutOffOutletV",0);
    //scalar weightOutletLV = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightOutletLV",0);

    /*forAll(costFunctionPatches,cI)
    {
        Foam::label patchI = mesh.boundaryMesh().findPatchID(costFunctionPatches[cI]);
        const Foam::fvPatch& currPatch = mesh.boundary()[patchI];

	forAll(currPatch,facei){
	    Foam::label faceCelli = currPatch.faceCells()[facei];
	    if (U[faceCelli].component(1) < cutOffOutletV){
		J += weightOutletLV*mag(U[faceCelli].component(1));
	    }
	}
    }*/

    forAll(costFunctionPatches,cI)
    {
        Foam::label patchI = mesh.boundaryMesh().findPatchID(costFunctionPatches[cI]);
        //const Foam::fvPatch& currPatch = mesh.boundary()[patchI];

	    /*forAll (currPatch,facei){
	        Foam::label faceCelli = currPatch.faceCells()[facei];
	        J -= (p.internalField()[faceCelli] + 0.5*rho.internalField()[faceCelli]*magSqr(U.internalField()[faceCelli]));
	    }*/

	    //J -= gSum(phi.boundaryField()[patchI]*(p.boundaryField()[patchI] + 0.5*magSqr(phi.boundaryField()[patchI]/currPatch.magSf())));
	    //J -= gSum(phi.boundaryField()[patchI]*(p.boundaryField()[patchI] + 0.5*rho.boundaryField()[patchI]*magSqr(U.boundaryField()[patchI])));

	    J -= gSum(
		    phi.boundaryField()[patchI]*p.boundaryField()[patchI]*pow(
		        (1+((gamma-1)/(2*gamma))*psi.boundaryField()[patchI]*magSqr(U.boundaryField()[patchI])),(gamma/(gamma-1))
		    )
	    );

	    /*outletTotalPressure += gSum(
		    p.boundaryField()[patchI]*pow(
		        (1+((gamma-1)/(2*gamma))*psi.boundaryField()[patchI]*magSqr(U.boundaryField()[patchI])),(gamma/(gamma-1))
		    )
	    );*/

	    //J -= gSum(p.boundaryField()[patchI]*(1+0.5*psi.boundaryField()[patchI]*magSqr(U.boundaryField()[patchI])));

	    /*forAll(currPatch,facei){
	        scalar MafI = mag(U.boundaryField()[patchI][facei])/cSound.boundaryField()[patchI][facei];
	        scalar psiI = psi.boundaryField()[patchI][facei];
	        scalar phiI = phi.boundaryField()[patchI][facei];
	        scalar pI = p.boundaryField()[patchI][facei];
	        scalar magSqrUI = magSqr(U.boundaryField()[patchI][facei]);
	        scalar rhoI = rho.boundaryField()[patchI][facei];

	        if (MafI >= 1){
		        J -= phiI*pI*pow((1 + (gamma-1)/(2*gamma)*psiI*magSqrUI),(gamma/(gamma-1)));
	        }else if (0.7 <= MafI && MafI < 1){
		        J -= phiI*pI*(1 + 0.5*psiI*magSqrUI);
	        }else{
		        J -= phiI*(pI + 0.5*rhoI*magSqrUI);
	        }
	    }*/
    }

    //Info << "\noutlet total pressure : " << outletTotalPressure << nl << endl;

    return J;
}

Foam::scalar CostFunction::eval_vorticity(const Foam::wordList &costFunctionPatches)
{
    const Foam::volVectorField& U = mesh.lookupObject<Foam::volVectorField>("U");

    scalar J = 0;
    
    volScalarField magVort = mag(fvc::curl(U));
    
    J = gSum(magVort);

    return J;
}

Foam::scalar CostFunction::eval_reacSonicTotalPressure(const Foam::wordList &costFunctionPatches)
{
    const Foam::volScalarField& p = mesh.lookupObject<Foam::volScalarField>("p");
    const Foam::surfaceScalarField& phi = mesh.lookupObject<Foam::surfaceScalarField>("phi");
    const Foam::volScalarField& psi = mesh.lookupObject<Foam::volScalarField>("thermo:psi");
    //const Foam::volScalarField& Cp = mesh.lookupObject<Foam::volScalarField>("thermo:cp");
    //const Foam::volScalarField& Cv = mesh.lookupObject<Foam::volScalarField>("thermo:cv");
    const Foam::volScalarField& gamma = mesh.lookupObject<Foam::volScalarField>("gamma");
    const Foam::volVectorField& U = mesh.lookupObject<Foam::volVectorField>("U");
    const Foam::volScalarField& cSound = mesh.lookupObject<Foam::volScalarField>("cSound");
    const Foam::volScalarField& rho = mesh.lookupObject<Foam::volScalarField>("rho");

    scalar J = 0;
    //scalar outletTotalPressure = 0;

    //const Foam::volScalarField& Ma = mag(U)/cSound;

    forAll(costFunctionPatches,cI)
    {
        Foam::label patchI = mesh.boundaryMesh().findPatchID(costFunctionPatches[cI]);
        const Foam::fvPatch& currPatch = mesh.boundary()[patchI];

	    /*J -= gSum(
		    p.boundaryField()[patchI]*pow(
		        (1+((gamma.boundaryField()[patchI]-1)/(2*gamma.boundaryField()[patchI]))*psi.boundaryField()[patchI]*magSqr(U.boundaryField()[patchI])),
		        (gamma.boundaryField()[patchI]/(gamma.boundaryField()[patchI]-1))
		    )
	        );*/

	    /*outletTotalPressure += gSum(
		    p.boundaryField()[patchI]*pow(
		        (1+((gamma.boundaryField()[patchI]-1)/(2*gamma.boundaryField()[patchI]))*psi.boundaryField()[patchI]*magSqr(U.boundaryField()[patchI])),
		        (gamma.boundaryField()[patchI]/(gamma.boundaryField()[patchI]-1))
		    )
	        );*/

	    forAll(currPatch,facei){
	        scalar MaI = mag(U.boundaryField()[patchI][facei])/cSound.boundaryField()[patchI][facei];
	        scalar psiI = psi.boundaryField()[patchI][facei];
	        scalar phiI = phi.boundaryField()[patchI][facei];
	        scalar pI = p.boundaryField()[patchI][facei];
	        scalar magSqrUI = magSqr(U.boundaryField()[patchI][facei]);
	        scalar rhoI = rho.boundaryField()[patchI][facei];
	        scalar gammaI = gamma.boundaryField()[patchI][facei];

	        if (MaI >= 1){
		        J -= phiI*pI*pow((1 + (gammaI-1)/(2*gammaI)*psiI*magSqrUI),(gammaI/(gammaI-1)));
	        }else if (0.7 <= MaI && MaI < 1){
		        J -= phiI*pI*(1 + 0.5*psiI*magSqrUI);
	        }else{
		        J -= phiI*(pI + 0.5*rhoI*magSqrUI);
	        }
	    }
    }

    //Info << "\noutlet total pressure : " << outletTotalPressure << nl << endl;

    return J;
}

Foam::scalar CostFunction::eval_sonicFlowPoint(const Foam::wordList &costFunctionPatches)
{
    const Foam::volScalarField& p = mesh.lookupObject<Foam::volScalarField>("p");
    const Foam::surfaceScalarField& phi = mesh.lookupObject<Foam::surfaceScalarField>("phi");
    const Foam::volScalarField& psi = mesh.lookupObject<Foam::volScalarField>("thermo:psi");
    const Foam::volScalarField& gamma = mesh.lookupObject<Foam::volScalarField>("gamma");
    const Foam::volVectorField& U = mesh.lookupObject<Foam::volVectorField>("U");
    const Foam::volScalarField& rho = mesh.lookupObject<Foam::volScalarField>("rho");
    const Foam::volScalarField& alpha = mesh.lookupObject<Foam::volScalarField>("alpha");
    const Foam::volScalarField& c = mesh.lookupObject<Foam::volScalarField>("c");

    scalar J = 0;
    //scalar totalPressureOut = 0;

    //scalar weightPoint = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightPoint",1); //CL
    //scalar weightOutlet = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightOutlet",1); //CL

    scalar targetMa =  mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("targetMa",1); //CL
    //scalar targetU =  mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("targetU",1); //CL

    scalar weightP = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightP",0); //CL
    //scalar weightGradP = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightGradP",0); //CL
    //scalar weightGradU = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightGradU",0); //CL
    //scalar weightDomainV = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightDomainV",0); //CL
    scalar weightAlpha = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightAlpha",0); //CL
    scalar weightTargetMa = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightTargetMa",0); //CL
    scalar weightOutletMa = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightOutletMa",0); //CL
    //scalar weightOutletMa = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightOutletMa",0); //CL
    //scalar weightOutletU = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightOutletU",0); //CL
    //scalar weightU = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("weightU",0); //CL

    scalar targetAlphaFrac = mesh.solutionDict().subDict("OPTIMIZATION").lookupOrDefault<scalar>("targetAlphaFrac",0.3);
    dimensionedScalar alphaMax(mesh.solutionDict().subDict("OPTIMIZATION").lookup("alphaMax")); //CL
    scalar alphaFrac = gSum(alpha)/(alphaMax.value()*alpha.size());

    //volTensorField gradU = Foam::fvc::grad(U);
    //volScalarField gradUx = gradU.component(0);
    //volScalarField gradUy = gradU.component(1);
    //volVectorField gradP = Foam::fvc::grad(p);

    /*const Foam::vector position = Foam::vector(
	    mesh.solutionDict().subDict("OPTIMIZATION").lookup("costFunctionPos")
	);*/

    //Foam::label pos = mesh.findCell(position);

    const Foam::vector p1 = Foam::vector(
            mesh.solutionDict().subDict("OPTIMIZATION").lookup("bottomLeft")
        );

    const Foam::vector p2 = Foam::vector(
            mesh.solutionDict().subDict("OPTIMIZATION").lookup("topRight")
        );

    //scalar sectionU = 0;
    //scalar sectionP = 0;
    scalar sectionMa = 0;
    //scalar outletMa = 0;
    //scalar outletU = 0;
    //scalar outletP = 0;
    scalar n = 0;
    //scalar domainGradP = 0;
    //scalar domainGradU = 0;
    //scalar domainV = 0;
    //scalar outletTotalPressure = 0;

    forAll(U,celli){
	    //domainGradP += (mag(gradP[celli].component(0)) + mag(gradP[celli].component(1)));
	    //domainGradU += (mag(gradUx[celli]) + mag(gradUy[celli]));
	    //domainV += mag(U[celli].component(1));

	    scalar xPos = mesh.C()[celli].component(0);
	    scalar yPos = mesh.C()[celli].component(1);

	    if (p1.component(0) <= xPos && xPos <= p2.component(0)){
	        if (p1.component(1) <= yPos && yPos <= p2.component(1)){
	            if (alpha[celli] == 0){
		            //sectionU += mag(U[celli]);
		            //sectionU += mag(U[celli].component(0));
		            //sectionP += p[celli];
		            //sectionMa += abs(U[celli].component(0))/cSound[celli];
		            sectionMa += mag(U[celli])/c[celli];
		            n += 1;
		        }else{
		            sectionMa += 0;
		            n += 1;
		        }
	        }
	    }
    }

    //scalar meanSectionU = sectionU/n;
    //scalar meanSectionP = sectionP/n;
    scalar meanSectionMa = sectionMa/n;
    //scalar outletMa = 0;

    /*forAll(costFunctionPatches,cI){
        Foam::label patchI = mesh.boundaryMesh().findPatchID(costFunctionPatches[cI]);
        const Foam::fvPatch& currPatch = mesh.boundary()[patchI];

	    forAll(currPatch,facei){
	        //outletMa += mag(U.boundaryField()[patchI][facei])/cSound.boundaryField()[patchI][facei];
	        //outletU += magSqr(U.boundaryField()[patchI][facei])*phi.boundaryField()[patchI][facei];
	        //outletP += gSum(p.boundaryField()[patchI])/currPatch.size();
	        
	        scalar MaI = mag(U.boundaryField()[patchI][facei])/cSound.boundaryField()[patchI][facei];
	        scalar psiI = psi.boundaryField()[patchI][facei];
	        scalar phiI = phi.boundaryField()[patchI][facei];
	        scalar pI = p.boundaryField()[patchI][facei];
	        scalar magSqrUI = magSqr(U.boundaryField()[patchI][facei]);
	        scalar rhoI = rho.boundaryField()[patchI][facei];
	        scalar gammaI = gamma.boundaryField()[patchI][facei];
	        
	        outletMa += MaI;

	        if (MaI >= 1){
		        outletTotalPressure += phiI*pI*pow((1 + (gammaI-1)/(2*gammaI)*psiI*magSqrUI),(gammaI/(gammaI-1)));
	        }else if (0.7 <= MaI && MaI < 1){
		        outletTotalPressure += phiI*pI*(1 + 0.5*psiI*magSqrUI);
	        }else{
		        outletTotalPressure += phiI*(pI + 0.5*rhoI*magSqrUI);
	        }
	    }
    }*/

    Info << "\nMean Section Mach Number : " << meanSectionMa << nl << endl;
    //Info << "Mean Section U : " << sectionU/n << nl << endl;
    Info << "Section Size : " << n << " Cells" << nl << endl;

    //J = -sectionU + weightP*sectionP + weightGradP*domainGradP + weightDomainV*domainV + weightAlphaFrac*pow(alphaFrac - targetAlphaFrac,2);
    //J = -weightU*sectionU + weightMa*pow(sectionMa/n - targetMa,2) - weightOutletU*outletU + weightAlpha*pow(alphaFrac - targetAlphaFrac,2);
    //J = weightTargetMa*pow(meanSectionMa - targetMa,2) + weightAlpha*pow(alphaFrac - targetAlphaFrac,2) - weightOutletMa*outletMa;
    J = weightTargetMa*pow(meanSectionMa - targetMa,2) + weightAlpha*pow(alphaFrac - targetAlphaFrac,2);
    //J = weightU*pow(sectionU/n - targetU,2) + weightP*sectionP;

    return J;
}
