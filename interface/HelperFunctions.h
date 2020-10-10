#ifndef HELPERFUNCTIONS_H
#define HELPERFUNCTIONS_H

#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TVectorD.h>

#include "Analysis/SelectorTools/interface/GoodParticle.h"

double JetSphericity(std::vector<GoodPart>& jets) {
    TMatrixDSym sphere(3);
    double pSum = 0;
    for(auto jet: jets) {
	for(int i = 1; i <= 3; i++) {
	    sphere(i-1,i-1) += pow(jet[i], 2);
	    pSum += pow(jet[i], 2);
	    for(int j = i+1; j <= 3; j++) {
		sphere(i-1,j-1) += jet[i]*jet[j];
		sphere(j-1,i-1) += jet[i]*jet[j];
	    }
	}
    }
    sphere *= 1/pSum;
    TMatrixDSymEigen eigenGetter(sphere);
    TVectorD eigen = eigenGetter.GetEigenValues();
    return 3./2*(eigen(1) + eigen(2));
}

double JetCentrality(std::vector<GoodPart>& jets, double HT) {
    double eTot = 0;
    for(auto jet: jets) {
	eTot += jet.E();
    }
    return HT/eTot;
}

#endif /* HELPERFUNCTIONS_H */

