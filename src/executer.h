
#ifndef SRC_EXECUTER_H_
#define SRC_EXECUTER_H_
#include "../extern/xtensor/include/xtensor/xarray.hpp"
#include "../extern/xtensor/include/xtensor/xadapt.hpp"

namespace openbps {

namespace executer {

void run_solver();

void init_solver();

void exponental(xt::xarray<double>& matrix, xt::xarray<double>& y);
 
void iterative(xt::xarray<double>& matrix, xt::xarray<double>& y);
void cram(xt::xarray<double>& matrix, xt::xarray<double>& y);
}

}





#endif /* SRC_EXECUTER_H_ */
