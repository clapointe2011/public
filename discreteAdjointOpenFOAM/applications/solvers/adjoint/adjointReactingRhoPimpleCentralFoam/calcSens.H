#include "fvCFD.H"

if (buildFrom && nIter == 1){

    sens == 0*sens;

    const fvPatchList& patches = mesh.boundary();
    forAll(patches, patchi){
	const fvPatch& currPatch = patches[patchi];
	if (isType<wallFvPatch>(currPatch)){
	    forAll(currPatch,facei){
		label faceCelli = currPatch.faceCells()[facei];
		sens[faceCelli] = dco::derivative(alpha[faceCelli]);
	    }
	}
    }
}/*else if (proximityDist && nIter > 1){
    forAll(sens,i){
	scalar xPosi = mesh.C()[i].component(0);
	scalar yPosi = mesh.C()[i].component(1);

	scalar dist = 1e6;
	scalar minDist = 1e6;

	forAll(alpha,j){
	    if (alpha[j] > 0 && i != j){
		scalar xPosj = mesh.C()[j].component(0);
		scalar yPosj = mesh.C()[j].component(1);
		dist = sqrt(pow(xPosj-xPosi,2) + pow(yPosj-yPosi,2));
	    }

	    if (dist < minDist){
		minDist = dist;
	    }
	}

	if (minDist < maxSeparation){
	    sens[i] = dco::derivative(alpha[i]);
	}else{
	    sens[i] = 0;
	}
    }
}*/else if (proximityNeighbor && nIter > 1){
    forAll(mesh.C(),celli){
	const labelList& neighborCells = mesh.cellCells()[celli];

	sens[celli] = 0;

	forAll(neighborCells,nI){
	    if (alpha[neighborCells[nI]] > 0){
		sens[celli] = dco::derivative(alpha[celli]);
		break;
	    }
	}
    }
}else{
    forAll(alpha,i){
	sens[i] = dco::derivative(alpha[i]);
    }
}

if (scaleByVol){
    #include "scaleByVol.H"
}
