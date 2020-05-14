#include "executer.h"

#include "chain.h"
#include "parse.h"
#include "configure.h"
#include "functionals.h"
#include "reactions.h"
#include "../extern/xtensor/include/xtensor/xarray.hpp"
#include "../extern/xtensor/include/xtensor/xadapt.hpp"

namespace openbps {

namespace executer {

void run_solver() {
//	std::vector<double> ytarget {collapsing(xval, yval, xtarget)};
	std::cout << "We star execution\n";
   



}

void exponental(xt::xarray<double>& matrix, xt::xarray<double>& y) {
           double dt {configure::timestep/configure::numstep};
			   std::cout << "MATRIX IS:" << std::endl;
			   matrix = xt::exp2(matrix * dt);
			   std::cout << "INIT OF CALCULATION IS: " <<  y.size() << std::endl;
			   for (int j = 0; j < y.size(); j++) {
			   		std::cout << chain.nuclides[j].name << " = " << y[j] << std::endl;

			   }
               for (int k = 0; k < configure::numstep; k++) {
            	   for (int j = 0; j < y.size(); j++) {
            		   for (int k = 0; k < y.size(); k++) {
            			    y(j) = y(j) * matrix(j, k);
            		   }
            	   }

			   }
}
 
void iterative(xt::xarray<double>& matrix, xt::xarray<double>& y){
           double dt {configure::timestep/configure::numstep};
			   std::vector<std::size_t> shape = { y.size() };
           xt::xarray<double> sigp = xt::zeros<double>(shape);
           xt::xarray<double> ro = xt::zeros<double>(shape);
           xt::xarray<double> roo = xt::zeros<double>(shape);
           xt::xarray<double> rrr = xt::zeros<double>(shape);
           xt::xarray<double> arr = xt::zeros<double>(shape);
           xt::xarray<double, 2> arrtemp = xt::zeros<double>(shape);
           xt::xarray<double> disr = xt::zeros<double>(shape);
           xt::xarray<double> rest = xt::zeros<double>(shape);
           xt::xarray<double> et = xt::ones<double>(shape);
           xt::xarray<double> er = xt::ones<double>(shape);
           size_t N {y.size()};
           double aa{0.0};
           for (size_t i=0; i < N; i++) {
	             sigp(i) = xt::sum(xt::view(matrix, 0, i));
                for (size_t j=0; j < N; j++) {
	                  matrix(i,j) = matrix(i,j)/sigp(i);
                }
           }
           arr = sigp * dt;
           disr = xt::exp2(-arr);
           rest = 1 - disr;
           //from fortran where (arr /= 0.0d-10) er = rest / arr  
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
                    for (int iparent=0; iparent < N; iparent++) {
	                     for (int ichild=0; ichild < N; ichild++) {
	                         arrtemp(i, j) = matrix(iparent, ichild)*rrr(iparent);
                        }
                     }
                     for (int ichild=0; ichild < N; ichild++) {
                         aa = xt::sum(xt::view(arrtemp, 1, ichild));
                         ro(ichild) += aa;
                     }
                     for (int iparent=0; iparent < N; iparent++) {
                         arr(iparent) = ro(iparent) - y(iparent);
                         rrr(iparent) = arr(iparent) * et(iparent);
                         y(iparent) += arr(iparent) * er(iparent);
                     }
                     ro = y;
                     //arr=0;
                     for (int iparent=0; iparent < N; iparent++) {
                         if (roo(iparent) > 0.0) {
                             if ( abs(1.0 - ro(iparent)/roo(iparent)) > epb) {
    proceed=true;
    break;
} else { proceed=false;}
                        }
                     }
                     roo = y;
               }
           }
}
void cram(xt::xarray<double>& matrix, xt::xarray<double>& y){
}

void init_solver() {

	//std::string input_file{"./Xmls/chain.xml"};
	pugi::xml_node chain_node{read_chain_xml(configure::chain_file)};

	Chain chain(chain_node);
	std::cout << "Simple reaction\n";
	    auto react_map = chain.form_reaction();
	    for(int i = 0; i < react_map["(n,gamma)"].size(); i++)
	    {
	        std::cout << react_map["(n,gamma)"][i].first << "  " << react_map["(n,gamma)"][i].second << "  "<< 1 << "  " <<std::endl;
	    }

	    std::cout << "NFY reaction\n";
	    auto map = chain.form_yield_map();
	    std::cout << "Energies of reaction\n";
	    for(auto iter = map.begin() ; iter != map.end() ; iter++)
	    {
	        std::cout << "Energy is " << iter->first << std::endl;
	    }
	    for(int i = 0; i < map[0.0253].size(); i++)
	    {
	        std::cout << map[0.0253][i][0] << "  " << map[0.0253][i][1] << "  "<< map[0.0253][i][2] << "  " <<std::endl;
	    }
	    read_reactions_xml();
       //forming from reactions INP
       std::vector<Materials> v = read_materials_from_reactions();
       form_materials_xml(v, XML_INP_PATH);
       //reading v from INP and form to OUT
       v = read_materials_from_inp(XML_INP_PATH);
       matchcompositions(v, compositions);
		for (auto& mat : v) {
			if (mat.name != "all") {
			   xt::xarray<double> mainarr;
			   std::vector<std::size_t> shape = { chain.nuclides.size() };
			   xt::xarray<double> y = make_concentration(chain, mat.namenuclides, mat.conc);
			   xt::xarray<double> y_old = xt::zeros<double>(shape);
			   y_old = y;
			   mainarr = form_matrix(chain, *mat.compos);
           case configure::mode {
               exponent: exponental(mainarr, y);
               iteration: iterative(mainarr, y);
               chebyshev: cram(mainarr, y);
               
           }
			   std::cout << "RESULT OF CALCULATION IS: " <<  y.size() << std::endl;
			   for (int j = 0; j < y.size(); j++) {
			   		std::cout << chain.nuclides[j].name << " = " << y[j] << std::endl;
			   }
			}

		}

}

}

}
#include "executer.h"

#include "chain.h"
#include "parse.h"
#include "configure.h"
#include "functionals.h"
#include "reactions.h"
#include "../extern/xtensor/include/xtensor/xarray.hpp"
#include "../extern/xtensor/include/xtensor/xadapt.hpp"

namespace openbps {

namespace executer {

void run_solver() {

	std::vector<double> xval = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
	std::vector<double> yval = {10.0, 22.0, 32.0, 33.0, 14.0, 85.0};
	std::vector<double> xtarget = {0.0, 6.0, 8.2};

//	std::vector<double> ytarget {collapsing(xval, yval, xtarget)};
	std::cout << "We star execution\n";



}

void init_solver() {

	//std::string input_file{"./Xmls/chain.xml"};
	pugi::xml_node chain_node{read_chain_xml(configure::chain_file)};

	Chain chain(chain_node);
	std::cout << "Simple reaction\n";
	    auto react_map = chain.form_reaction();
	    for(int i = 0; i < react_map["(n,gamma)"].size(); i++)
	    {
	        std::cout << react_map["(n,gamma)"][i].first << "  " << react_map["(n,gamma)"][i].second << "  "<< 1 << "  " <<std::endl;
	    }

	    std::cout << "NFY reaction\n";
	    auto map = chain.form_yield_map();
	    std::cout << "Energies of reaction\n";
	    for(auto iter = map.begin() ; iter != map.end() ; iter++)
	    {
	        std::cout << "Energy is " << iter->first << std::endl;
	    }
	    for(int i = 0; i < map[0.0253].size(); i++)
	    {
	        std::cout << map[0.0253][i][0] << "  " << map[0.0253][i][1] << "  "<< map[0.0253][i][2] << "  " <<std::endl;
	    }
	    read_reactions_xml();
		for (auto& compos : compositions) {
			if (compos.name != "all") {
			   xt::xarray<double> mainarr;
			   std::vector<std::size_t> shape = { chain.nuclides.size() };
			   xt::xarray<double> y = make_concentration(chain, compos.namenuclides, compos.conc);
			   xt::xarray<double> y_old = xt::zeros<double>(shape);
			   y_old = y;
			   double dt {configure::timestep/configure::numstep};
			   mainarr = form_matrix(chain, compos);
			   std::cout << "MATRIX IS:" << std::endl;
			   mainarr = xt::exp2(mainarr * dt);
			   std::cout << "INIT OF CALCULATION IS: " <<  y.size() << std::endl;
			   for (int j = 0; j < y.size(); j++) {
			   		std::cout << chain.nuclides[j].name << " = " << y[j] << std::endl;

			   }
               for (int k = 0; k < configure::numstep; k++) {
            	   for (int j = 0; j < y.size(); j++) {
            		   for (int k = 0; k < y.size(); k++) {
            			    y(j) = y(j) * mainarr(j, k);
            		   }
            	   }

			   }
			   std::cout << "RESULT OF CALCULATION IS: " <<  y.size() << std::endl;
			   for (int j = 0; j < y.size(); j++) {
			   		std::cout << chain.nuclides[j].name << " = " << y[j] << std::endl;
			   }
			}

		}

}

}

}
