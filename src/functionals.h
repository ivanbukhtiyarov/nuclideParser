#ifndef SRC_FUNCTIONALS_H_
#define SRC_FUNCTIONALS_H_

#include <vector>

namespace openbps {

std::vector<double> collapsing(const std::vector<double>& xval,
		                       const std::vector<double>& yval,
							   const std::vector<double>& xtarget);

std::vector<std::pair<int, double>>
translating(const std::vector<double>& xval, const std::vector<double>& xtarget);



std::vector<double> transition(const std::vector<double>& xval,
		                       const std::vector<double>& yval,
							   const std::vector<double>& xtarget);

} //namespace openbps
#endif /* SRC_FUNCTIONALS_H_ */
