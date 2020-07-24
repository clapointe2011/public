if (smoothResField && nIter % smoothFreq == 0){
    forAll(mesh.C(),celli){
	    const labelList& neighborCells = mesh.cellCells()[celli];
	    scalar sumLocalAlpha = 0;
	    bool isSurrounded = true;

	    forAll(neighborCells,nI){
	        sumLocalAlpha += alpha[neighborCells[nI]];

	        if (alpha[neighborCells[nI]] == 0){
		        isSurrounded = false;
	        }
	    }

	    scalar meanLocalAlpha = sumLocalAlpha/neighborCells.size();

	    if (isSurrounded){
	        if (fillGaps && alpha[celli] == 0){
		        smoothAlpha[celli] = meanLocalAlpha;
	        }else{
		        smoothAlpha[celli] = SAweight*alpha[celli] + (1-SAweight)*meanLocalAlpha;
	        }
	    }else if (!isSurrounded && meanLocalAlpha > SAfactor*alphaMax.value()){
	        smoothAlpha[celli] = SAweight*alpha[celli] + (1-SAweight)*meanLocalAlpha;
	    }else{
	        smoothAlpha[celli] = 0;
	    }
    }

    alpha == smoothAlpha;
}

if (limitMinAlpha){
    forAll(alpha,i){
	    if (alpha[i] < MAfactor*alphaMax.value()){
	        alpha[i] = 0;
	    }
    }
}
