#include <cmath>
#include "chain.h"
#include "functionals.h"
#include "parse.h"
#include "configure.h"
#include "reactions.h"
#include "../extern/xtensor/include/xtensor/xarray.hpp"
#include "../extern/xtensor/include/xtensor/xadapt.hpp"
//
namespace openbps {

  std::vector<std::pair<int, std::string>> Chain::form_idx_name() {
      std::vector<std::pair<int, std::string>> out;
      std::pair<int, std::string> elem;

      for (int i = 0; i < nuclides.size(); i++) {
          elem.first = i;
          elem.second = nuclides[i].name;
          out.push_back(elem);
        //  std::cout << "idx" << elem.first << "name" << elem.second << std::endl;
      }
      return out;
  }

  std::vector<std::pair<int, double>> Chain::form_idx_lambda() {
      std::vector<std::pair<int, double>> out;
      std::pair<int, double> elem;

      for (int i = 0; i < nuclides.size(); i++) {
          elem.first = i;
          elem.second = nuclides[i].half_life;
          out.push_back(elem);
          // std::cout << "idx" << elem.first << "Lamda " << elem.second << std::endl;
      }
      return out;
  }

  std::vector<std::vector<double>> Chain::form_idx_decay() {
      std::vector<std::vector<double>> out;
      std::vector<double> elem;

      for (int i = 0; i < nuclides.size(); i++) {
          for (int j = 0; j < nuclides[i].decay_arr.size(); j++) {
              elem.push_back(i);
              elem.push_back(name_idx[nuclides[i].decay_arr[j].target]);
              elem.push_back(nuclides[i].decay_arr[j].branching_ratio);
              out.push_back(elem);
              elem.clear();
          }
      }
      return (out);
  }
/*
* to test the following parser use:
*
*
  auto react_map = chain.form_reaction();
  for(int i = 0; i < react_map["(n,gamma)"].size(); i++)
  {
      std::cout << react_map["(n,gamma)"][i].first << "  " << react_map["(n,gamma)"][i].second << "  "<< 1 << "  " <<std::endl;
  }
  output:
0  1  1
1  2  1
7  8  1
17  18  1
18  19  1
28  29  1
40  41  1
41  42  1
...
*/
  std::map<std::string, std::vector<std::pair<int, int>>> Chain::form_reaction() {

      std::map<std::string, std::vector<std::pair<int, int>>> out;
      std::vector<std::pair<int, int>> n_2n;
      std::vector<std::pair<int, int>> n_3n;
      std::vector<std::pair<int, int>> n_a;
      std::vector<std::pair<int, int>> n_gamma;
      std::vector<std::pair<int, int>> n_p;
      std::pair<int, int> elem;

      out.insert({"(n,2n)", n_2n});
      out.insert({"(n,3n)", n_3n});
      out.insert({"(n,gamma)", n_gamma});
      out.insert({"(n,a)", n_a});
      out.insert({"(n,p)", n_p});

      for (int i = 0; i < nuclides.size(); i++) {
          for (int j = 0; j < nuclides[i].reaction_arr.size(); j++) {
              elem.first = i;
              elem.second = name_idx[nuclides[i].reaction_arr[j].target];
              out[nuclides[i].reaction_arr[j].type].push_back(elem);
          }
      }
      return (out);
  }

/*
* To test this parser (only Pu239 have 2000000.0):
  auto map = chain.form_yield_map();
  for(int i = 0; i < map[2000000.0].size(); i++)
  {
      std::cout << map[2000000.0][i][0] << "  " << map[2000000.0][i][1] << "  "<< map[2000000.0][i][2] << "  " <<std::endl;
  }
output
...
3557  1346  0
3557  1350  0
3557  1351  0
3557  1352  6.12916e-14
3557  1353  1.65977e-13
3557  1354  2.47967e-12
and so on
*/

  std::map<double, std::vector<std::vector<double>>> Chain::form_yield_map() {

      std::map<double, std::vector<std::vector<double>>> out;
      std::vector<std::vector<double>> insertion;
      std::vector<double> elem;

      for(int i = 0; i < nuclides.size(); i++) {

          for(int j = 0 ; j < nuclides[i].nfy.energies.size(); j++) {

              if(out.count(nuclides[i].nfy.energies[j]) == 0)
                  out.insert({nuclides[i].nfy.energies[j], insertion});
          }
      }
      for(int i = 0; i < nuclides.size(); i++) {

          for(int j = 0 ; j < nuclides[i].nfy.energies.size(); j++) {

              auto n_map = nuclides[i].nfy.yield_arr[j].product_data;
              for(auto& item : n_map)
              {
                  elem.push_back(i);
                  elem.push_back(name_idx[item.first]);
                  elem.push_back(item.second);
                  out[nuclides[i].nfy.yield_arr[j].energy].push_back(elem);
                  elem.clear();
              }
          }
      }
      return out;
  }

  bool comp (std::pair <double, double> a,std::pair <double, double> b) {
    return a.first < b.first;
  }

  std::pair<std::vector<double>, std::vector<double>>
	Chain::get_yield_map_(size_t father, const std::string& daughter) {
          std::pair<std::vector<double>, std::vector<double>> out;
          std::vector<std::pair<double, double>> fracenergy_;
          std::vector<double> energies;
          std::vector<double> fractions;

          for(int j = 0 ; j < nuclides[father].nfy.energies.size(); j++) {

              auto n_map = nuclides[father].nfy.yield_arr[j].product_data;
              auto it = n_map.find(daughter);
              if (it != n_map.end()) {
                  fracenergy_.push_back(
                  std::make_pair(nuclides[father].nfy.energies[j],
                      		   it->second));
             }
          }
          std::sort(fracenergy_.begin(), fracenergy_.end(), comp);
          for (auto it=fracenergy_.begin(); it != fracenergy_.end(); it++){
              energies.push_back(it->first);
              fractions.push_back(it->second);
          }
           out = std::make_pair(energies, fractions);

          return out;
 }

  xt::xarray<double> form_matrix(Chain& chainer, Materials& mat) {
  	size_t NN {chainer.nuclides.size()};
  	int k {0};
  	xt::xarray<int>::shape_type shape = {NN, NN};
  	xt::xarray<double> result = xt::zeros<double>(shape);
  	std::pair<std::vector<double>, std::vector<double>> pair1, pair2;
  	std::vector<std::string> fplist;
  	std::vector<double> weight;
  	double decay_ {0.0};
        if (mat.compos != nullptr) {
            pair2 = mat.compos->get_fluxenergy();
        }
  	for (int i = 0; i < NN; i++) {

          if (chainer.nuclides[i].half_life > 0) { //~0

          	decay_ = log(2.0) / chainer.nuclides[i].half_life;
          	result(i, i) = -decay_;
          	for (int j = 0; j < chainer.nuclides[i].decay_arr.size(); j++) {
          		if (chainer.nuclides[i].decay_arr[j].target != std::string("Nothing")) {
          			k = chainer.name_idx[chainer.nuclides[i].decay_arr[j].target];
          			result(k, i) += decay_* chainer.nuclides[i].decay_arr[j].branching_ratio;
          		}
          	}

          }
          if ((mat.compos != nullptr) && (mat.compos->energy_number > 0)) {// 0 --  4-5 secs

             for (auto& obj : mat.compos->xslib) {//0 - 4,5 secs per iteration
          	   if (obj.xsname == chainer.nuclides[i].name) { // average 1 sec max 4 s

                     double rr {0.0};
                     for (auto& r: obj.rxs) {//exactly ~0
                     			rr += r;

                     }
                     for (int j = 0; j < chainer.nuclides[i].reaction_arr.size(); j++) {// 0-3sec
                         if (obj.xstype == chainer.nuclides[i].reaction_arr[j].type) {//0-3sec avg 1 sec
                      	   if (obj.xstype != "fission") {
                      	       k = chainer.name_idx[chainer.nuclides[i].reaction_arr[j].target];
                      	       result(k, i) += rr * PWD;//*branching ratio
                      	   } else {// ELSE works 0.8-3 s

                      		   for (int j = 0 ; j < chainer.nuclides[i].nfy.energies.size(); j++)  { //each iter 0,3 - 1,8 s
                      		        auto n_map = chainer.nuclides[i].nfy.yield_arr[j].product_data;
                      		        for (auto& item : n_map) {//0,0007 - 0,006 sec
                      		            if (std::find(fplist.begin(), fplist.end(), item.first) == fplist.end()) {

                      		            	pair1 = chainer.get_yield_map_(i, item.first);// 1-5 ms
                      		            	k = chainer.name_idx[item.first];
                      		                fplist.push_back(item.first);
                      		                double br {0.0};
                      		                double norm {0.0};
                                              if (weight.empty())
                      		                    weight = transition(pair2.first, pair2.second, pair1.first);//~0
                      		                for (int l = 0; l < weight.size(); l++) {
                                                  br += weight[l] * pair1.second[l];
                                                  norm += weight[l];
                      		                }

                      		                result(k, i) = br / norm * rr * PWD;
                      		                norm = result(k, i);
                                                norm = 1.;
                      		            }
                      		        }
                                  weight.clear();
                      		   }
                      	   fplist.clear();
                      	   }
                         }
                     }
                     result(i, i) -= rr * PWD;
          	   }
             }

          }
  	}
  	return result;

  }
  //
  xt::xarray<double> form_sigp(Chain& chainer, Materials& mat) {
  	size_t NN {chainer.nuclides.size()};
  	int k {0};
  	xt::xarray<int>::shape_type shape = {NN};
  	xt::xarray<double> result = xt::zeros<double>(shape);
  	double decay_ {0.0};

  	for (int i = 0; i < NN; i++) {

          if (chainer.nuclides[i].half_life > 0) { //~0

          	decay_ = log(2.0) / chainer.nuclides[i].half_life;
          	result(i) = decay_;
          }
          if ((mat.compos != nullptr) && (mat.compos->energy_number > 0)) {// 0 --  4-5 secs

             for (auto& obj : mat.compos->xslib) {//0 - 4,5 secs per iteration
          	   if (obj.xsname == chainer.nuclides[i].name) { // average 1 sec max 4 s
                     double rr {0.0};
                     for (auto& r: obj.rxs) {rr += r;}
                     result(i) += rr * PWD;
                   }

             }
          }
        }

  	return result;

  }
  //
  xt::xarray<double> make_concentration(Chain& chainer, std::vector<std::string>& nameconc,
  		                              std::vector<double>& ro) {
  	size_t NN {chainer.nuclides.size()};
  	xt::xarray<int>::shape_type shape = {NN};
  	xt::xarray<double> result = xt::zeros<double>(shape);
  	//for (int i = 0; i < nameconc.size(); i++) {
        for (int i = 0; i < ro.size(); i++) {
  		if (chainer.name_idx.find(nameconc[i]) != chainer.name_idx.end()) {
  			result[chainer.name_idx[nameconc[i]]] = ro[i];
  		}

  	}
  	return result;

  }
 
  xt::xarray<double> form_dmatrix(Chain& chainer, Materials& mat) {
  	size_t NN {chainer.nuclides.size()};
  	int k {0};
  	xt::xarray<int>::shape_type shape = {NN, NN};
  	xt::xarray<double> result = xt::zeros<double>(shape);
  	std::pair<std::vector<double>, std::vector<double>> pair1, pair2;
  	std::vector<std::string> fplist;
  	std::vector<double> weight;
  	double decay_ {0.0};
        if (mat.compos != nullptr) {
            pair2 = mat.compos->get_fluxenergy();
        }
  	for (int i = 0; i < NN; i++) {

          if (chainer.nuclides[i].half_life > 0) { //~0

          	decay_ = log(2.0) / (chainer.nuclides[i].half_life * chainer.nuclides[i].half_life) * chainer.nuclides[i].d_half_life;//-
          	//decay_ = chainer.nuclides[i].d_half_life / chainer.nuclides[i].half_life;//+
          	result(i, i) = -decay_;
          	for (int j = 0; j < chainer.nuclides[i].decay_arr.size(); j++) {
          		if (chainer.nuclides[i].decay_arr[j].target != std::string("Nothing")) {
          			k = chainer.name_idx[chainer.nuclides[i].decay_arr[j].target];
          			result(k, i) += decay_* chainer.nuclides[i].decay_arr[j].branching_ratio + 
                                                chainer.nuclides[i].decay_arr[j].d_branching_ratio * log(2.0) / chainer.nuclides[i].half_life;
          		}
          	}

          }
          if ((mat.compos != nullptr) && (mat.compos->energy_number > 0)) {// 0 --  4-5 secs

             for (auto& obj : mat.compos->xslib) {//0 - 4,5 secs per iteration
          	   if (obj.xsname == chainer.nuclides[i].name) { // average 1 sec max 4 s

                     double rr {0.0};
                     double drr {0.0};
                     for (auto& r: obj.rxs) {//exactly ~0
                     			rr += r;

                     }
                     for (auto& r: obj.d_rxs) {//exactly ~0
                     			drr += r;

                     }
                     for (int j = 0; j < chainer.nuclides[i].reaction_arr.size(); j++) {// 0-3sec
                         if (obj.xstype == chainer.nuclides[i].reaction_arr[j].type) {//0-3sec avg 1 sec
                      	   if (obj.xstype != "fission") {
                      	       k = chainer.name_idx[chainer.nuclides[i].reaction_arr[j].target];
                      	       result(k, i) += drr * PWD;//*branching ratio
                      	   } else {// ELSE works 0.8-3 s
                                   if (obj.xsname =="U235") std::cout <<"UU235 drr = " << drr <<std::endl;
                      		   for (int j = 0 ; j < chainer.nuclides[i].nfy.energies.size(); j++)  { //each iter 0,3 - 1,8 s
                      		        auto n_map = chainer.nuclides[i].nfy.yield_arr[j].product_data;
                      		        for (auto& item : n_map) {//0,0007 - 0,006 sec
                      		            if (std::find(fplist.begin(), fplist.end(), item.first) == fplist.end()) {

                      		            	pair1 = chainer.get_yield_map_(i, item.first);// 1-5 ms
                      		            	k = chainer.name_idx[item.first];
                      		                fplist.push_back(item.first);
                      		                double br {0.0};
                      		                double norm {0.0};
                                              if (weight.empty())
                      		                    weight = transition(pair2.first, pair2.second, pair1.first);//~0
                      		                for (int l = 0; l < weight.size(); l++) {
                                                  br += weight[l] * pair1.second[l];
                                                  norm += weight[l];
                      		                }

                      		                result(k, i) = br / norm * drr * PWD;
                      		                norm = result(k, i);
                                                norm = 1.;
                      		            }
                      		        }
                                  weight.clear();
                      		   }
                      	   fplist.clear();
                      	   }
                         }
                     }
                     result(i, i) -= drr * PWD;
          	   }
             }

          }
  	}
  	return result;
  }
  xt::xarray<double> form_dsigp(Chain& chainer, Materials& mat) {
        size_t NN {chainer.nuclides.size()};
  	int k {0};
  	xt::xarray<int>::shape_type shape = {NN};
  	xt::xarray<double> result = xt::zeros<double>(shape);
  	double decay_ {0.0};

  	for (int i = 0; i < NN; i++) {

          if (chainer.nuclides[i].half_life > 0) { //~0

          	
          	decay_ = log(2.0) / (chainer.nuclides[i].half_life * chainer.nuclides[i].half_life) * chainer.nuclides[i].d_half_life;//-
          	//decay_ = chainer.nuclides[i].d_half_life / chainer.nuclides[i].half_life;//+
          	result(i) = decay_;
          }
          if ((mat.compos != nullptr) && (mat.compos->energy_number > 0)) {// 0 --  4-5 secs

             for (auto& obj : mat.compos->xslib) {//0 - 4,5 secs per iteration
          	   if (obj.xsname == chainer.nuclides[i].name) { // average 1 sec max 4 s
                     double drr {0.0};
                     for (auto& r: obj.d_rxs) {drr += r;}
                     result(i) += drr * PWD;
                   }

             }
          }
        }

  	return result;

  }

  // Read a chain from xml
  pugi::xml_node read_chain_xml(const std::string& filename) {
      //pugi::xml_document doc;
      auto result = configure::docx.load_file(filename.c_str());
      if (!result) {
          std::cerr << "Error: file not found!" << std::endl;
      }

      pugi::xml_node chain_node = configure::docx.child("depletion_chain");
      return chain_node;
  }


} // namespace close
