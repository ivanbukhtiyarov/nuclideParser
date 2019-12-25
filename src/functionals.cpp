#include "functionals.h"
#include <iostream>
#include <vector>

namespace openbps {

std::vector<double> collapsing(const std::vector<double>& xval,
		               const std::vector<double>& yval,
                               const std::vector<double>& xtarget){
	std::cout << "Collapse\n";

	if ((xval.size() < 2) || (xtarget.size() < 2)) {
		std::cerr << "Incorrect vector sizes to collapse it\n";
	}

	if (yval.size() != xval.size() - 1){
		std::cerr << "Incorrect vector sizes to collapse it\n";
	}

	std::vector<double> ytarget(xtarget.size() + 1);
	if (xtarget[xtarget.size()-1] < xval[0]) {
		for (int i = 0; i < xtarget.size(); i++){
			ytarget[i] = yval[0];
		}
		return ytarget;
	}

	if (xtarget[0] > xval[xval.size() - 1]) {
		for (int i = 0; i < xtarget.size(); i++){
			ytarget[i] = yval[0];
		}
		return ytarget;
	}
    int ic {0};
    int jc {1};
    double val {yval[0]};
    double normx {0.0};
    double bx {xtarget[0]};
    for (int i = 1; i < xval.size(); i++){
    	while ((jc < xtarget.size()) && (xtarget[jc] <= xval[i])) {
    		normx += val * (xtarget[jc] - bx);
    		ytarget[jc-1] = normx / (xtarget[jc] - xtarget[jc-1]);
    		bx = xtarget[jc];
    		normx = 0.0;
    		jc += 1;
    	}
    	if (jc >= xtarget.size()) {
    		return ytarget;
    	} else {
    		normx += val * (xval[i] - bx);
    	    bx = xval[i];
    	    if (i < xval.size() - 1)  val = yval[i];
    	}

    }
    if (jc < xtarget.size()) {
        for (int j = jc; j < xtarget.size(); j++) {
        	normx += val * (xtarget[j] - bx);
            ytarget[j-1] = normx / (xtarget[j] - xtarget[j-1]);
            bx = xtarget[j];
            normx = 0.0;
        }
    }
    return ytarget;

}


std::vector<std::pair<int, double>>
translating(const std::vector<double>& xval, const std::vector<double>& xtarget) {
	std::cout << "Translate\n";
	std::vector<std::pair<int, double>> ytranslater;
	std::pair <int, double> tpair (0, 2.30);
	ytranslater.push_back(tpair);
	return ytranslater;
}

}




