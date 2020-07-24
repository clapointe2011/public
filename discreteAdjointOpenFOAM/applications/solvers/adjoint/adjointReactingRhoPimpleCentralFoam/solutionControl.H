//check for convergence

Info << "\nchecking for optimization convergence" << nl << endl;
Info << "cost function value : " << J << nl << endl;
Info << "cost function delta : " << abs(J - prevJ) << nl << endl;
Info << "alpha total : " << gSum(alpha) << nl << endl;
Info << "alpha delta : " << abs(gSum(alpha) - prevAlpha) << nl << endl;

delta = abs(gSum(alpha)-prevAlpha);
tol = abs(relTol*gSum(alpha));

//delta = abs(J-prevJ) + abs(gSum(alpha) - prevAlpha);
//tol = abs(relTol*(J+gSum(alpha)));

//delta = abs(J-prevJ);
//tol = abs(relTol*prevJ);

if (delta<tol && nIter > 1){
    /*if (perturbSolution && pIter <= perturbCutoff){

	Info << "cost function delta within tolerance" << nl << endl;
	Info << "perturbing solution" << nl << endl;

	forAll(mesh.C(),celli){

	    //CL zero perturbed velocity field
	    perturbUField[celli] = vector::zero;

	    Random rand(123456); 

	    const labelList& neighborCells = mesh.cellCells()[celli];
	    scalar localUx = 0;
	    scalar localUy = 0;

	    forAll(neighborCells,nI){
		localUx += U[neighborCells[nI]].component(0);
		localUy += U[neighborCells[nI]].component(1);
	    }

	    localUx += U[celli].component(0);
	    localUy += U[celli].component(1);

	    scalar meanLocalUx = localUx/(neighborCells.size()+1);
	    scalar meanLocalUy = localUy/(neighborCells.size()+1);

	    //CL new perturbed velocity field value is sum of original velocity component and a random amound of a local mean
	    perturbUField[celli].component(0) = U[celli].component(0) + PSfactor*rand.GaussNormal()*meanLocalUx;
	    perturbUField[celli].component(1) = U[celli].component(1) + PSfactor*rand.GaussNormal()*meanLocalUy;
	}

	forAll(U,i){
	    U[i] = perturbUField[i];
	}

	//CL write perturbed field
	perturbUField.write();

	//CL store old values
	prevJ = CostFunction(mesh).eval();
	prevAlpha = gSum(alpha);

	Info << "solution perturbed" << nl << endl;
	Info << "continuing optimization" << nl << endl;

    }else{
	#include "postProcAlpha.H"

	Info << "cost function delta within tolerance" << nl << endl;
	Info << "stopping optimization" << nl << endl;

	if (passiveFinalization){
	    dco::ga1s<double>::global_tape->switch_to_passive();
	    #include "passiveFinalization.H"
	}

	//CL write final fields
	U.write(); //CL
	p.write(); //CL
	T.write(); //CL
	alpha.write(); //CL

	break;
    }*/

    #include "postProcAlpha.H"

    Info << "alpha delta within tolerance" << nl << endl;
    Info << "stopping optimization" << nl << endl;

    if (passiveFinalization){
	    dco::ga1s<double>::global_tape->switch_to_passive();
	    #include "passiveFinalization.H"
    }

    //CL write final fields
    runTime.setTime(0,(label)0); //CL
    runTime.writeNow(); //CL

    CostFunction(mesh).eval();

    break;
}else if (nIter == optEndIter){
    #include "postProcAlpha.H"

    Info << "max optimization iteration reached" << nl << endl;
    Info << "stopping optimization" << nl << endl;

    if (passiveFinalization){
	    dco::ga1s<double>::global_tape->switch_to_passive();
	    #include "passiveFinalization.H"
    }

    //CL write final fields
    runTime.setTime(0,(label)0); //CL
    runTime.writeNow(); //CL

    CostFunction(mesh).eval();

}else{
    //CL store old values
    prevJ = CostFunction(mesh).eval();
    prevAlpha = gSum(alpha);

    Info << "alpha delta not within tolerance" << nl << endl;
    Info << "continuing optimization" << nl << endl;
}
