#include "reactions.h"
#include <vector>
#include <map>
#include <iostream>
#include <string.h>
#include "configure.h"
#include "functionals.h"
#include "parse.h"
#include "../extern/pugiData/pugixml.h"
#include "../extern/xtensor/include/xtensor/xarray.hpp"
#include "../extern/xtensor/include/xtensor/xadapt.hpp"

namespace openbps {

std::vector<Composition> compositions;
std::map<std::string, int> composmap;
size_t indexall;

Composition::Composition(pugi::xml_node node) {

	if (check_for_node(node, "name")) {
	    this->name = node.attribute("name").value();
	}
	if (check_for_node(node, "energy")) {
		for (pugi::xml_node tool : node.children("energy")) {
            size_t cng = std::stoi(tool.attribute("ng").value());
            this->energies_.insert({cng, get_node_array<double>(tool, "energy")});
            this->energy_number++;
		}
	}
	if (check_for_node(node, "flux")) {
	    for (pugi::xml_node tool : node.children("flux")) {
	        this->flux_= get_node_array<double>(tool, "flux");
		}
	}
        if (check_for_node(node, "dflux")) {
	    for (pugi::xml_node tool : node.children("dflux")) {
	        this->d_flux_= get_node_array<double>(tool, "dflux");
		}
	}

	if (check_for_node(node, "spectrum")) {
	    for (pugi::xml_node tool : node.children("spectrum")) {
	        size_t cng = std::stoi(tool.attribute("ng").value());
	        this->spectrum_= get_node_array<double>(tool, "spectrum");
		}
	}

        if (check_for_node(node, "dspectrum")) {
	    for (pugi::xml_node tool : node.children("dspectrum")) {
	        this->d_spectrum_= get_node_array<double>(tool, "dspectrum");
		}
	}

	if (check_for_node(node, "namenuclides")) {
	    this->namenuclides = get_node_array<std::string>(node, "namenuclides");
	}

	if (check_for_node(node, "conc")) {
		this->conc = get_node_array<double>(node, "conc");
	}

	if (check_for_node(node, "xslibs")) {
           std::string rxs {node.child("xslibs").attribute("typex").value()};
           for (pugi::xml_node tool : node.child("xslibs").children("xslib")) {
		this->xslib.push_back(parse_xs_xml_(tool, rxs, 0));
           }
           for (pugi::xml_node tool : node.child("xslibs").children("dxslib")) {
		this->xslib.push_back(parse_xs_xml_(tool, rxs, 1));
           }
	}
        
	this->nuclide_number = this->conc.size();
	this->spectrum_.size() > 0 ?
	    this->energy_number = this->spectrum_.size() :
		this->energy_number = this->flux_.size();
}


s_xs_ Composition::parse_xs_xml_(pugi::xml_node node, std::string& rxs, int mode) {
    s_xs_ result;    
    result.xstype = node.attribute("reaction").value();
    result.xsname = node.attribute("name").value();
    if (rxs == "cs") {
        mode == 0 ? 
                   result.xs_ = get_node_array<double>(node, "xslib") : 
                   result.d_xs_ = get_node_array<double>(node, "dxslib");
    } else {
        mode == 0 ? 
    	           result.rxs = get_node_array<double>(node, "xslib"):
                   result.d_rxs = get_node_array<double>(node, "dxslib");
    }
    return result;
}

void Composition::depcopymap_(std::map<size_t, std::vector<double>>& fmap,
		        std::map<size_t, std::vector<double>>& smap) {

	std::map<size_t, std::vector<double>>::iterator it;
	if (!smap.empty()) {
		   for ( it= smap.begin(); it != smap.end(); ++it ) {
			   if (fmap.find(it->first) == fmap.end()){
				   fmap[it->first] = it->second;
			   }

		   }
	   }

}

void Composition::deploy_all(Composition& externcompos) {

   if (externcompos.name == this->name) {
	   return;
   } else {

	   this->depcopymap_(this->energies_, externcompos.energies_);
	   if (!externcompos.spectrum_.empty() && this->spectrum_.empty()) {
	       this->spectrum_.resize(externcompos.spectrum_.size());
	       std::copy(externcompos.spectrum_.begin(),externcompos.spectrum_.end(),
	   	   	 this->spectrum_.begin());
	   }
	   if (!externcompos.flux_.empty() && this->flux_.empty()) {
	   		   this->flux_.resize(externcompos.flux_.size());
	   		   std::copy(externcompos.flux_.begin(),externcompos.flux_.end(),
	   				     this->flux_.begin());
	   }
           if (!externcompos.d_spectrum_.empty() && this->d_spectrum_.empty()) {
	       this->d_spectrum_.resize(externcompos.d_spectrum_.size());
	       std::copy(externcompos.d_spectrum_.begin(),externcompos.d_spectrum_.end(),
	   	   	 this->d_spectrum_.begin());
	   }
	   if (!externcompos.d_flux_.empty() && this->d_flux_.empty()) {
	   		   this->d_flux_.resize(externcompos.d_flux_.size());
	   		   std::copy(externcompos.d_flux_.begin(),externcompos.d_flux_.end(),
	   				     this->d_flux_.begin());
	   }
	   if (!externcompos.namenuclides.empty() && this->namenuclides.empty()) {
		   this->namenuclides.resize(externcompos.namenuclides.size());
		   std::copy(externcompos.namenuclides.begin(),externcompos.namenuclides.end(),
				     this->namenuclides.begin());
	   }
	   if (!externcompos.conc.empty() && this->conc.empty()) {
		   this->conc.resize(externcompos.conc.size());
		   std::copy(externcompos.conc.begin(),externcompos.conc.end(),
					 this->conc.begin());
	   }
	   if (!externcompos.xslib.empty() && this->xslib.empty()){
		   this->xslib.resize(externcompos.xslib.size());
		   std::copy(externcompos.xslib.begin(),externcompos.xslib.end(),
		   			 this->xslib.begin());

	   }
   }
}

void  Composition::get_reaction() {

	for (auto ixs = this->xslib.begin(); ixs != this->xslib.end(); ixs++) {
	if (ixs->rxs.empty()){
		ixs->rxs.resize(1);
		ixs->rxs[0] = 0.0;
        if (!ixs->xs_.empty()) {
        	int ng = ixs->xs_.size();
        	std::vector<double> cuflux(ng, 1.0);
        	if (this->flux_.size() > 0) {
        	if (this->flux_.size() == ng) {
        		cuflux = this->flux_;
        	} else {
        		cuflux = collapsing(this->energies_[this->flux_.size()],
        		                    this->flux_,
					    this->energies_[ng]);
        	}
        	}
            for (int i = 0; i != ng; i++) {
            	ixs->rxs[0] += cuflux[i] * ixs->xs_[i];
            }

        }
	}
        if (ixs->d_rxs.empty()){
		ixs->d_rxs.resize(1);
		ixs->d_rxs[0] = 0.0;
        if (!ixs->d_xs_.empty()) {
        	int ng = ixs->d_xs_.size();
        	std::vector<double> cuflux(ng, 1.0);
        	std::vector<double> d_cuflux(ng, 0.0);
        	if (this->flux_.size() > 0) {
        	if (this->flux_.size() == ng) {
        		cuflux = this->flux_;
        	} else {
        		cuflux = collapsing(this->energies_[this->flux_.size()],
        		                    this->flux_,
					    this->energies_[ng]);
        	}
        	}
        	int ng_d = ixs->xs_.size();
        	if (this->d_flux_.size() > 0) {
                if (this->d_flux_.size() == ng_d) {
        		d_cuflux = this->d_flux_;
        	} else {
        		d_cuflux = collapsing(this->energies_[this->d_flux_.size()],
        		                    this->d_flux_,
					    this->energies_[ng_d]);
        	}
        	}
            for (int i = 0; i != ng; i++) {
            	ixs->d_rxs[0] += cuflux[i] * ixs->d_xs_[i];
            }
            for (int i = 0; i != ng_d; i++) {
                ixs->d_rxs[0] += d_cuflux[i] * ixs->xs_[i];
            }

        }
	}

}
}

bool compfe (std::pair <double, double> a,std::pair <double, double> b) {
  return a.first < b.first;
}

std::pair<std::vector<double>, std::vector<double>> Composition::get_fluxenergy() {
	std::pair<std::vector<double>, std::vector<double>> result;
	size_t ng;
	double fluxnorm {0.0};
	std::vector<std::pair<double, double>> forsort_;
	if (spectrum_.empty()) {
		std::for_each(flux_.begin(), flux_.end(), [&] (double n) {
			fluxnorm += n;
		});
		for (auto& f: flux_) {
			spectrum_.push_back(f/fluxnorm);
		}
	}
	ng = spectrum_.size();

    result = std::make_pair(energies_[ng], spectrum_);
    return result;
}

void read_reactions_xml() {
	pugi::xml_document doc;
	auto result = doc.load_file(configure::reaction_file.c_str());
	if (!result) {
	    std::cout << "Warning: file reactions.xml not found!" << std::endl;
            return;
	}
	pugi::xml_node root_node = doc.child("compositions");

	std::cout << "I' m in reactions.xml parser" << std::endl;

	for (pugi::xml_node tool : root_node.children("composit")) {
		compositions.push_back(Composition(tool));
		int index = compositions.size() - 1;
		if (compositions[index].name == "all") {
			indexall = index;
		}
		composmap.insert({compositions[index].name, index});
	}

	for (std::vector<Composition>::iterator it = compositions.begin();
	     it != compositions.end() ;++it) {
        it->deploy_all(compositions[indexall]);
        it->get_reaction();
	}


}

}


