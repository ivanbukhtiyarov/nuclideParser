#include "executer.h"
#include <cmath>
#include "chain.h"
#include "parse.h"
#include "configure.h"
#include "functionals.h"
#include "reactions.h"
#include <iostream>
#include <fstream>
#include <complex>
#include "../extern/xtensor/include/xtensor/xarray.hpp"
#include "../extern/xtensor/include/xtensor/xadapt.hpp"
#include "../extern/xtensor/include/xtensor/xmath.hpp"
#include "../extern/xtensor/include/xtensor/xview.hpp"
#include "../extern/xtensor/include/xtensor/xbuilder.hpp"
#include "../extern/xtensor/include/xtensor/xcomplex.hpp"
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <vector>


typedef Eigen::SparseMatrix<double> SpMat;

using namespace std::complex_literals;
std::string XML_INP_PATH = "./Xmls/materials.xml";

namespace openbps {

namespace executer {

//cram
xt::xarray<std::complex<double>> theta16 {
    +3.509103608414918 + 8.436198985884374i,
    +5.948152268951177 + 3.587457362018322i,
    -5.264971343442647 + 16.22022147316793i,
    +1.419375897185666 + 10.92536348449672i,
    +6.416177699099435 + 1.194122393370139i,
    +4.993174737717997 + 5.996881713603942i,
    -1.413928462488886 + 13.49772569889275i,
    -10.84391707869699 + 19.27744616718165i};

xt::xarray<std::complex<double>> alpha16{
    +5.464930576870210e+3 - 3.797983575308356e+4i,
    +9.045112476907548e+1 - 1.115537522430261e+3i,
    +2.344818070467641e+2 - 4.228020157070496e+2i,
    +9.453304067358312e+1 - 2.951294291446048e+2i,
    +7.283792954673409e+2 - 1.205646080220011e+5i,
    +3.648229059594851e+1 - 1.155509621409682e+2i,
    +2.547321630156819e+1 - 2.639500283021502e+1i,
    +2.394538338734709e+1 - 5.650522971778156e+0i};

xt::xarray<double> theta_48r{
    -4.465731934165702e+1, -5.284616241568964e+0,
    -8.867715667624458e+0, +3.493013124279215e+0,
    +1.564102508858634e+1, +1.742097597385893e+1,
    -2.834466755180654e+1, +1.661569367939544e+1,
    +8.011836167974721e+0, -2.056267541998229e+0,
    +1.449208170441839e+1, +1.853807176907916e+1,
    +9.932562704505182e+0, -2.244223871767187e+1,
    +8.590014121680897e-1, -1.286192925744479e+1,
    +1.164596909542055e+1, +1.806076684783089e+1,
    +5.870672154659249e+0, -3.542938819659747e+1,
    +1.901323489060250e+1, +1.885508331552577e+1,
    -1.734689708174982e+1, +1.316284237125190e+1};

xt::xarray<double> theta_48i{
    +6.233225190695437e+1, +4.057499381311059e+1,
    +4.325515754166724e+1, +3.281615453173585e+1,
    +1.558061616372237e+1, +1.076629305714420e+1,
    +5.492841024648724e+1, +1.316994930024688e+1,
    +2.780232111309410e+1, +3.794824788914354e+1,
    +1.799988210051809e+1, +5.974332563100539e+0,
    +2.532823409972962e+1, +5.179633600312162e+1,
    +3.536456194294350e+1, +4.600304902833652e+1,
    +2.287153304140217e+1, +8.368200580099821e+0,
    +3.029700159040121e+1, +5.834381701800013e+1,
    +1.194282058271408e+0, +3.583428564427879e+0,
    +4.883941101108207e+1, +2.042951874827759e+1};

xt::xarray<double> alpha_48r{
    +6.387380733878774e+2, +1.909896179065730e+2,
    +4.236195226571914e+2, +4.645770595258726e+2,
    +7.765163276752433e+2, +1.907115136768522e+3,
    +2.909892685603256e+3, +1.944772206620450e+2,
    +1.382799786972332e+5, +5.628442079602433e+3,
    +2.151681283794220e+2, +1.324720240514420e+3,
    +1.617548476343347e+4, +1.112729040439685e+2,
    +1.074624783191125e+2, +8.835727765158191e+1,
    +9.354078136054179e+1, +9.418142823531573e+1,
    +1.040012390717851e+2, +6.861882624343235e+1,
    +8.766654491283722e+1, +1.056007619389650e+2,
    +7.738987569039419e+1, +1.041366366475571e+2};

xt::xarray<double> alpha_48i{
    -6.743912502859256e+2, -3.973203432721332e+2,
    -2.041233768918671e+3, -1.652917287299683e+3,
    -1.783617639907328e+4, -5.887068595142284e+4,
    -9.953255345514560e+3, -1.427131226068449e+3,
    -3.256885197214938e+6, -2.924284515884309e+4,
    -1.121774011188224e+3, -6.370088443140973e+4,
    -1.008798413156542e+6, -8.837109731680418e+1,
    -1.457246116408180e+2, -6.388286188419360e+1,
    -2.195424319460237e+2, -6.719055740098035e+2,
    -1.693747595553868e+2, -1.177598523430493e+1,
    -4.596464999363902e+3, -1.738294585524067e+3,
    -4.311715386228984e+1, -2.777743732451969e+2};
xt::xarray<std::complex<double>> theta48{theta_48r+theta_48i*1.i};
xt::xarray<std::complex<double>> alpha48 {alpha_48r+alpha_48i*1.i};;
double alpha480{2.258038182743983e-47};
double alpha160{2.124853710495224e-16};
//cram

void run_solver() {
    std::cout << "We star execution\n";
}

void exponental(xt::xarray<double>& matrix, xt::xarray<double>& y) {
    double dt {configure::timestep/configure::numstep};
    std::vector<std::vector<double>> out;
    if (configure::outwrite) {
	configure::dumpoutput.clear();
	configure::dumpoutput.resize(configure::numstep);
        for (int k = 0; k < configure::numstep; k++) {
            configure::dumpoutput[k].resize(y.size());
        }
    }
    xt::xarray<double> matrix2, ro;
    ro = y * 0.0;
    matrix2 = xt::exp(matrix * dt);
    for (int k = 0; k < configure::numstep; k++) {
        for (int j = 0; j < y.size(); j++) {

             for (int i = 0; i < y.size(); i++) {
                  if (i==j) {
                      ro(j) += y(i) * matrix2(j, i);
                  } else {
                      ro(j) += y(i) *  (matrix2(j, i) - 1);
                  }
             }
             if ((std::isnan(ro(j))) || ( std::isinf(ro(j)))) ro(j) = 0.0;

        }
        y = ro;
        ro = ro * 0.0;
        if (configure::outwrite) {
            for (int i = 0; i < y.size(); i++) {
        	configure::dumpoutput[k][i] = y(i);
            }
        }

    }
}

void iterative(xt::xarray<double>& matrix, xt::xarray<double>& sigp, xt::xarray<double>& y){
           double dt {configure::timestep/configure::numstep};
	   std::vector<std::size_t> shape = { y.size() };
	   std::vector<std::size_t> dshape = { y.size(), y.size() };
           //xt::xarray<double> sigp = xt::zeros<double>(shape);
           xt::xarray<double> ro = xt::zeros<double>(shape);
           xt::xarray<double> roo = xt::zeros<double>(shape);
           xt::xarray<double> rrr = xt::zeros<double>(shape);
           xt::xarray<double> arr = xt::zeros<double>(shape);
           xt::xarray<double> arrtemp = xt::zeros<double>(dshape);//?!
           xt::xarray<double> disr = xt::zeros<double>(shape);
           xt::xarray<double> rest = xt::zeros<double>(shape);
           xt::xarray<double> et = xt::ones<double>(shape);
           xt::xarray<double> er = xt::ones<double>(shape);
           size_t N {y.size()};
           double aa {0.0};
           if (configure::outwrite) {
               configure::dumpoutput.clear();
               configure::dumpoutput.resize(configure::numstep);
               for (int k = 0; k < configure::numstep; k++) {
        	   configure::dumpoutput[k].resize(y.size());
               }
           }
           for (size_t i=0; i < N; i++) {
               matrix(i, i) = 0.0;
	       //sigp(i) = xt::sum(xt::col(matrix, i))(0);
           }
           arr = sigp * dt;
           disr = xt::exp(-arr);
           rest = 1 - disr;
           for (size_t i=0; i < N; i++) {
               if (rest(i) < 1.e-10) rest(i) = arr(i);
               if (arr(i) > 0.0) er(i) = rest(i) / arr(i);
               for (size_t j=0; j < N; j++) {
               	    if (sigp(j) > 0) matrix(i,j) = matrix(i,j) / sigp(j);

               }
           }
           et = 1.0 - er;
           roo = y;
           for (int k = 0; k < configure::numstep; k++) {
               bool proceed {true};
               int t = 0;
               rrr = y * rest;       //!< rest
               y = y * disr;         //!< disappearance
               ro = y;
               while (proceed) {
	            t++;
                    for (size_t iparent=0; iparent < N; iparent++) {
	                     for (size_t ichild=0; ichild < N; ichild++) {
	                         arrtemp(ichild, iparent) = matrix(ichild, iparent) * rrr(iparent);
                        }
                     }
                     for (size_t ichild=0; ichild < N; ichild++) {
                         aa = xt::sum(xt::row(arrtemp, ichild))(0);
                         ro(ichild) += aa;
                     }
                     arr = ro - y;
                     rrr = arr * et;
                     y += arr * er;
                     ro = y;

                     for (int iparent=0; iparent < N; iparent++) {
                         if (roo(iparent) > 0.0) {
                             if ( abs(1.0 - ro(iparent)/roo(iparent)) > configure::epb) {
                                  proceed=true;
                                  break;
                             } else {
                        	 proceed=false;
                             }
                        }
                     }
                     if (t > 25) {

                         std::cout << "Exceed :number of iteration " << std::endl;
                         proceed=false;
                     }
                     roo = y;
               }
               if (configure::outwrite) {
                   for (int i = 0; i < y.size(); i++) {
                       configure::dumpoutput[k][i] = y(i);
                   }
              }
           }//numstep
}
//
void diterative(xt::xarray<double>& matrix, xt::xarray<double>& sigp, xt::xarray<double>& y,
                xt::xarray<double>& dmatrix, xt::xarray<double>& dsigp, xt::xarray<double>& dy) {
           double dt {configure::timestep/configure::numstep};
	   std::vector<std::size_t> shape = { y.size() };
	   std::vector<std::size_t> dshape = { y.size(), y.size() };
           //xt::xarray<double> sigp = xt::zeros<double>(shape);
           xt::xarray<double> ro = xt::zeros<double>(shape);
           xt::xarray<double> dro = xt::zeros<double>(shape);
           xt::xarray<double> roo = xt::zeros<double>(shape);
           xt::xarray<double> rrr = xt::zeros<double>(shape);
           xt::xarray<double> drr = xt::zeros<double>(shape);
           xt::xarray<double> arr = xt::zeros<double>(shape);
           xt::xarray<double> arrtemp = xt::zeros<double>(dshape);
           xt::xarray<double> drrtemp = xt::zeros<double>(dshape);
           xt::xarray<double> disr = xt::zeros<double>(shape);
           xt::xarray<double> rest = xt::zeros<double>(shape);
           xt::xarray<double> et = xt::ones<double>(shape);
           xt::xarray<double> er = xt::ones<double>(shape);
           xt::xarray<double> ds = xt::ones<double>(shape);
           xt::xarray<double> dr = xt::ones<double>(shape);
           size_t N {y.size()};
           double aa {0.0};
           if (configure::outwrite) {
               configure::dumpoutput.clear();
               configure::dumpoutput.resize(configure::numstep);
               for (int k = 0; k < configure::numstep; k++) {
        	   configure::dumpoutput[k].resize(y.size());
               }
           }
           for (size_t i=0; i < N; i++) {
               matrix(i, i) = 0.0;
               dmatrix(i, i) = 0.0;
	       //sigp(i) = xt::sum(xt::col(matrix, i))(0);
           }
           arr = sigp * dt;
           disr = xt::exp(-arr);
           rest = 1 - disr;
           for (size_t i=0; i < N; i++) {
               if (rest(i) < 1.e-10) rest(i) = arr(i);
               if (arr(i) > 0.0) er(i) = rest(i) / arr(i);
               for (size_t j=0; j < N; j++) {
               	    if (sigp(j) > 0) matrix(i,j) = matrix(i,j) / sigp(j);
                    if (sigp(j) > 0) dmatrix(i,j) = dmatrix(i,j) / sigp(j);

               }
           }
           et = 1.0 - er;
           roo = y;
           for (int k = 0; k < configure::numstep; k++) {
               bool proceed {true};
               int t = 0;
               rrr = y * rest;       //!< rest
               drr = dy * rest;
               dro = disr * dy  + dsigp * dt;
               /*
               for (size_t iparent=0; iparent < N; iparent++) {
                   std :: cout << "DRO: " << dro(iparent) << std::endl;
                   std :: cout << "DY: " << dy(iparent) << std::endl;
                   std :: cout << "DSIGP: " << dsigp(iparent) << std::endl;
               }
               */
               dy = dro;
               y = y * disr;         //!< disappearance
               ro = y;
               while (proceed) {
	            t++;
                    for (size_t iparent=0; iparent < N; iparent++) {
	                     for (size_t ichild=0; ichild < N; ichild++) {
	                         arrtemp(ichild, iparent) = matrix(ichild, iparent) * rrr(iparent);
                                 drrtemp(ichild, iparent) = matrix(ichild, iparent) * drr(iparent) +  
                                                            dmatrix(ichild, iparent) * rrr(iparent);
                        }
                     }
                     for (size_t ichild=0; ichild < N; ichild++) {
                         aa = xt::sum(xt::row(arrtemp, ichild))(0);
                         ro(ichild) += aa;
                         dro(ichild) += xt::sum(xt::row(drrtemp, ichild))(0);
                     }
                     arr = ro - y;
                     ds = dro - dy;
                     rrr = arr * et;
                     drr = ds * et;
                     y += arr * er;
                     dy += ds * er;
                     ro = y;
                     dro = dy;
                     for (int iparent=0; iparent < N; iparent++) {
                         if (roo(iparent) > 0.0) {
                             if ( abs(1.0 - ro(iparent)/roo(iparent)) > configure::epb) {
                                  proceed=true;
                                  break;
                             } else {
                        	 proceed=false;
                             }
                        }
                     }
                     if (t > 25) {

                         std::cout << "Exceed :number of iteration " << std::endl;
                         proceed=false;
                     }
                     roo = y;
               }
               if (configure::outwrite) {
                   for (int i = 0; i < y.size(); i++) {
                       configure::dumpoutput[k][i] = y(i);
                   }
              }
           }//numstep
}

//
void diterative2(xt::xarray<double>& matrix, xt::xarray<double>& sigp, xt::xarray<double>& y,
                xt::xarray<double>& dmatrix, xt::xarray<double>& dsigp, xt::xarray<double>& dy) {
           double dt {configure::timestep/configure::numstep};
	   std::vector<std::size_t> shape = { y.size() };
	   std::vector<std::size_t> dshape = { y.size(), y.size() };
           //xt::xarray<double> sigp = xt::zeros<double>(shape);
           xt::xarray<double> ro = xt::zeros<double>(shape);
           xt::xarray<double> dro = xt::zeros<double>(shape);
           xt::xarray<double> roo = xt::zeros<double>(shape);
           xt::xarray<double> rrr = xt::zeros<double>(shape);
           xt::xarray<double> drr = xt::zeros<double>(shape);
           xt::xarray<double> arr = xt::zeros<double>(shape);
           xt::xarray<double> arrtemp = xt::zeros<double>(dshape);
           xt::xarray<double> drrtemp = xt::zeros<double>(dshape);
           xt::xarray<double> disr = xt::zeros<double>(shape);
           xt::xarray<double> rest = xt::zeros<double>(shape);
           xt::xarray<double> et = xt::ones<double>(shape);
           xt::xarray<double> er = xt::ones<double>(shape);
           xt::xarray<double> es = xt::ones<double>(shape);
           xt::xarray<double> ds = xt::ones<double>(shape);
           xt::xarray<double> dr = xt::ones<double>(shape);
           xt::xarray<double> ddr = xt::ones<double>(shape);
           size_t N {y.size()};
           double aa {0.0};
           double dd {0.0};
           if (configure::outwrite) {
               configure::dumpoutput.clear();
               configure::dumpoutput.resize(configure::numstep);
               for (int k = 0; k < configure::numstep; k++) {
        	   configure::dumpoutput[k].resize(y.size());
               }
           }
           for (size_t i=0; i < N; i++) {
               matrix(i, i) = 0.0;
               dmatrix(i, i) = 0.0;
	       //sigp(i) = xt::sum(xt::col(matrix, i))(0);
           }
           arr = sigp * dt;
           disr = xt::exp(-arr);
           rest = 1 - disr;
           for (size_t i=0; i < N; i++) {
               if (rest(i) < 1.e-10) rest(i) = arr(i);
               if (arr(i) > 0.0) er(i) = rest(i) / arr(i);
               for (size_t j=0; j < N; j++) {
               	    if (sigp(j) > 0) matrix(i,j) = matrix(i,j) / sigp(j);
                    if (sigp(j) > 0) dmatrix(i,j) = dmatrix(i,j) / sigp(j);

               }
           }
           et = 1.0 - er;
           es = (er + disr) * dsigp;
           
           roo = y;
           for (int k = 0; k < configure::numstep; k++) {
               bool proceed {true};
               int t = 0;
               rrr = y * rest;       //!< rest
               drr = dy * rest;
               ddr = (sigp * dt * dsigp + rest * dy);
               dr = y * disr * ddr;
               dy = dy + sigp * dt * dsigp;
               ds = y * disr * dy;
               /*
               for (size_t iparent=0; iparent < N; iparent++) {
                   std :: cout << "DRO: " << dro(iparent) << std::endl;
                   std :: cout << "DY: " << dy(iparent) << std::endl;
                   std :: cout << "DSIGP: " << dsigp(iparent) << std::endl;
               }
               */
               dro = dy;
               y = y * disr;         //!< disappearance
               ro = y;
               while (proceed) {
	            t++;
                    for (size_t iparent=0; iparent < N; iparent++) {
	                     for (size_t ichild=0; ichild < N; ichild++) {
	                         arrtemp(ichild, iparent) = matrix(ichild, iparent) * rrr(iparent);
                                 if (t == 0) {
                                     drrtemp(ichild, iparent) = dmatrix(ichild, iparent) + dy(iparent);
                                 } else {
                                     drrtemp(ichild, iparent) = dmatrix(ichild, iparent) + dr(iparent);
                                     
                                 }
                                 drrtemp(ichild, iparent) = arrtemp(ichild, iparent) * drrtemp(ichild, iparent);
                        }
                     }
                     for (size_t ichild=0; ichild < N; ichild++) {
                         aa = xt::sum(xt::row(arrtemp, ichild))(0);
                         dd = xt::sum(xt::row(drrtemp, ichild))(0);
                         ro(ichild) += aa;
                         dr(ichild) = aa * es(ichild);
                         ds(ichild) += dd * er(ichild) + dr(ichild);
                         dr(ichild) += dd * et(ichild); 
                     }
                     arr = ro - y;
                     
                     rrr = arr * et;
                     
                     y += arr * er;
                     
                     
                     for (int iparent=0; iparent < N; iparent++) {
                         if (y(iparent) < 0.0) y(iparent) = 0.0;
                         if (arr(iparent) > 0.0) dr(iparent) = dr(iparent) / arr(iparent);
                         if (dr(iparent) < 0.0) dr(iparent) = 0.0;
                         if (y(iparent) > 0.0) dy(iparent) = ds(iparent) ;// y(iparent);
                         //if (y(iparent) > 0.0) dy(iparent) = ds(iparent) / y(iparent);
                     }
                     ro = y;
                     for (int iparent=0; iparent < N; iparent++) {
                         if (roo(iparent) > 0.0) {
                             if ( abs(1.0 - ro(iparent)/roo(iparent)) > configure::epb) {
                                  proceed=true;
                                  break;
                             } else {
                        	 proceed=false;
                             }
                        }
                     }
                     if (t > 25) {

                         std::cout << "Exceed :number of iteration " << std::endl;
                         proceed=false;
                     }
                     roo = y;
               }
               if (configure::outwrite) {
                   for (int i = 0; i < y.size(); i++) {
                       configure::dumpoutput[k][i] = y(i);
                   }
              }
           }//numstep
}

void cram(xt::xarray<double>& matrix, xt::xarray<double>& y,
                                      xt::xarray<std::complex<double>>& alpha,
                                      xt::xarray<std::complex<double>>& theta,
                                      int order, double alpha0) {

    double dt {configure::timestep/configure::numstep};
    size_t n {y.size()};
    if (configure::outwrite) {
	configure::dumpoutput.clear();
	configure::dumpoutput.resize(configure::numstep);
        for (int k = 0; k < configure::numstep; k++) {
            configure::dumpoutput[k].resize(y.size());
        }
    }
    std::vector<std::size_t> shape = {y.size()};
    xt::xarray<std::complex<double>> zy(shape);
    xt::xarray<double> ident = xt::eye(y.size());
    matrix = matrix * dt;
    Eigen::VectorXd ytemp(n);
    Eigen::VectorXcd x(n), b(n), zytemp(n);
    Eigen::SparseMatrix<std::complex<double>> A(n,n);
    std::complex<double> c(1,0);
    A.reserve(n * 10);
    Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>, Eigen::COLAMDOrdering<int> >   solver;
    // fill A and b;
    for (size_t i=0; i < n; i++) {
         ytemp(i) = y(i);
         zytemp(i) = c * y(i);
         for (size_t j=0; j < n; j++) {
             if (matrix(i, j) != 0.) A.insert(i, j) = c * matrix(i, j);
         }
    }
    A.makeCompressed();

    for (int t=0; t < configure::numstep; t++) {

        for (int it=0; it < order; it++) {
            std :: cout << "Number of iterantions is "<<it<< std::endl;
            for (size_t j=0; j < n; j++) {
                 A.coeffRef(j,j) = c * matrix(j, j) - theta(it);
            }
            // Compute the ordering permutation vector from the structural pattern of A
            solver.analyzePattern(A);
            // Compute the numerical factorization
            solver.factorize(A);
            //Use the factors to solve the linear system
            zytemp.real() = ytemp;
            x = solver.solve(zytemp);

            x = alpha(it) * x;
	    ytemp += 2 * x.real();
            //y += 2 * real(alpha(it) * xt::linalg::solve(matrix - theta(it) * ident, zy));
        }

        ytemp = ytemp * alpha0;
        for (size_t i=0; i < n; i++)
            y(i) = ytemp(i);
        if (configure::outwrite) {
            for (int i = 0; i < y.size(); i++) {
                 configure::dumpoutput[t][i] = y(i);
            }
        }
    }

}

void init_solver() {

	pugi::xml_node chain_node{read_chain_xml(configure::chain_file)};

	Chain chain(chain_node);
	    read_reactions_xml();
            //forming from reactions INP
            std::vector<Materials> v = read_materials_from_inp(configure::inmaterials_file);
            matchcompositions(compositions, v);
            
	    for (auto& mat : v) {
	         if (mat.name != "all") {
		     xt::xarray<double> mainarr, dmainarr;
		     std::vector<std::size_t> shape = { chain.nuclides.size() };
		     xt::xarray<double> y = make_concentration(chain, mat.namenuclides, mat.conc);
                     if (configure::uncertantie_mod) std :: cout << "UNCERTATNTIES MOD!!!!!" << std :: endl;
                     xt::xarray<double> dy = make_concentration(chain, mat.namenuclides, mat.d_conc);
                     if (configure::uncertantie_mod) std :: cout << "UNCERTATNTIES MOD!!!!!" << std :: endl;
		     xt::xarray<double> y_old = xt::zeros<double>(shape);
                     xt::xarray<double> dy_old = xt::zeros<double>(shape);
		     y_old = y;
                     dy_old = dy;
	             mainarr = form_matrix(chain, mat);
                     switch (configure::calcmode) {
                       case configure::Mode::exponent:
                	 exponental(mainarr, y);
		         break;
                       case configure::Mode::iteration:
                         {
                          if (configure::uncertantie_mod) std :: cout << "1. out: out uncertanties mode" << std :: endl;
                          xt::xarray<double> sigp {form_sigp(chain, mat)} ;
                          if (configure::uncertantie_mod) {
                              xt::xarray<double> dsigp {form_dsigp(chain, mat)};
                              dmainarr = form_dmatrix(chain, mat);
                              if (configure::uncertantie_mod) std :: cout << "2. out: out uncertanties mode" << std :: endl;
                              diterative2(mainarr, sigp, y,
                                        dmainarr, dsigp, dy);
                          } else {
                	      iterative(mainarr, sigp, y);
                          }
                         }
		         break;
                       case configure::Mode::chebyshev:
                         if (configure::order == 8) {
                	     cram(mainarr, y, alpha16, theta16, configure::order, alpha160);
                         } else {
                             cram(mainarr, y, alpha48, theta48, configure::order, alpha480);
                         }
                	 break;
                     }
		     
	             for (size_t j = 0; j < y.size(); j++) {
			 if (y[j] > 0 && dy[j]/y[j] > 1.0) std::cout << chain.nuclides[j].name << " = " << y[j] <<" with error = " <<dy[j]<< std::endl;
			 if (configure::rewrite) mat.add_nuclide(chain.nuclides[j].name, y[j]);
			 if (configure::rewrite && configure::uncertantie_mod) mat.add_nuclide(chain.nuclides[j].name, dy[j], true);
		     }
                     std::cout << "RESULT OF CALCULATION IS: " <<  dy.size() << std::endl;
                     std::cout << "RESULT OF CALCULATION IS: " <<  y.size() << std::endl;
	             if (configure::outwrite) {
                                std::cout << "CHECK IN: " <<  y.size() << std::endl;
	             		std::ofstream myfile("outlog.csv", std::ofstream::app);
	             		if (myfile.is_open()) {
	             		    myfile << mat.name;
	             		    myfile << "\n";
                                    myfile << "Number" << ";" << "Act"<<";"<< "Q, Mev"<<"\n";
	             		    for (size_t t = 0; t < configure::numstep; t++) {
                                         double actval {0.0};
	             		         double qval {0.0};
                                         double acterr {0.0};
	             		         double qerr {0.0};
	             		         for (size_t j = 0; j < y.size(); j++) {
                                             //std :: cout << "FIND RES "<<chain.nuclides[j].name.rfind("U")<<std::endl;
	             		              if ((chain.nuclides[j].half_life > 0)) {
	             		        	  actval += log(2.0) / chain.nuclides[j].half_life * configure::dumpoutput[t][j] * 1.e+24;
	             		        	  qval += log(2.0) / chain.nuclides[j].half_life * configure::dumpoutput[t][j] *
                                                                                                   chain.nuclides[j].decay_energy * 1.e+24 *
                                                                                                   1.e-6;
                                                  if (configure::uncertantie_mod && y[j] > 0) {
                                                      /*acterr += log(2.0) / (chain.nuclides[j].half_life * chain.nuclides[j].half_life) *
                                                                chain.nuclides[j].d_half_life * configure::dumpoutput[t][j] * 1.e+24;*/
                                                      acterr += log(2.0) / chain.nuclides[j].half_life * dy[j] * 1.e+24;
                                                      qerr   += (log(2.0) / (chain.nuclides[j].half_life * chain.nuclides[j].half_life) *
                                                                chain.nuclides[j].d_half_life * configure::dumpoutput[t][j] * 1.e+24 +
                                                                log(2.0) / chain.nuclides[j].half_life * dy[j] * 1.e+24) * chain.nuclides[j].decay_energy * 1.e-6;
                                                  } 
                                                  //qval += configure::dumpoutput[t][j];
	             		              }
	             		          }
                                          if (!configure::uncertantie_mod) {
	             		              myfile << t << ";" << actval<<";" << qval << "\n";
                                          } else {
                                              myfile << t << ";" << actval<<";" << qval<< ";" << acterr<<";" << qerr << "\n";
                                          } 
	             		    }
	             		    myfile.close();
                                    std::cout << "CHECK OUT: " <<  y.size() << std::endl;
	             		}
	             	    }
		}
	    }

	    if (configure::rewrite) {
		form_materials_xml(v, configure::outmaterials_file);
	    }
            


}

}//namespace executer

}//namespace openbps
