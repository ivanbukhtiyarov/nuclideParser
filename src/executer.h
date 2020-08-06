#ifndef SRC_EXECUTER_H_
#define SRC_EXECUTER_H_
#include <complex>
#include "../extern/xtensor/include/xtensor/xarray.hpp"
#include "../extern/xtensor/include/xtensor/xadapt.hpp"
#include "../extern/xtensor/include/xtensor/xcomplex.hpp"
#include "../extern/xtensor/include/xtensor/xbuilder.hpp"
//cram

namespace openbps {

namespace executer {
//cram
extern xt::xarray<std::complex<double>> theta16;
extern xt::xarray<std::complex<double>> alpha16;
extern xt::xarray<std::complex<double>> theta48;
extern xt::xarray<std::complex<double>> alpha48;
extern xt::xarray<double> theta_48r;
extern xt::xarray<double> theta_48i;
extern xt::xarray<double> alpha_48r;
extern xt::xarray<double> alpha_48i;
extern double alpha480;
extern double alpha160;
//cram
void run_solver();

void init_solver();

void exponental(xt::xarray<double>& matrix, xt::xarray<double>& y);

void iterative(xt::xarray<double>& matrix, xt::xarray<double>& sigp, xt::xarray<double>& y);

void cram(xt::xarray<double>& matrix, xt::xarray<double>& y, 
                                      xt::xarray<std::complex<double>>& alpha, 
                                      xt::xarray<std::complex<double>>& theta, 
                                      int order, double alpha0);
}

}




#endif /* SRC_EXECUTER_H_ */
