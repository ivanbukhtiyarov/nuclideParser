#include "reactions.h"
#include <vector>
#include <map>
#include <iostream>
#include <string.h>
#include "configure.h"
#include "nuclide_class.h"
#include "../extern/pugiData/pugixml.h"
#include "../extern/xtensor/xarray.hpp"
#include "../extern/xtensor/xadapt.hpp"

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
	        size_t cng = std::stoi(tool.attribute("ng").value());
	        this->flux_.insert({cng, get_node_array<double>(tool, "flux")});
		}
	}
	if (check_for_node(node, "spectrum")) {
	    for (pugi::xml_node tool : node.children("spectrum")) {
	        size_t cng = std::stoi(tool.attribute("ng").value());
	        this->spectrum_.insert({cng, get_node_array<double>(tool, "spectrum")});
		}
	}

	if (check_for_node(node, "namenuclides")) {
	    this->namenuclides_ = get_node_array<std::string>(node, "namenuclides");
	}

	if (check_for_node(node, "conc")) {
		this->conc_ = get_node_array<double>(node, "conc");
	}

	if (check_for_node(node, "xslibs")) {
		this->xslib_.push_back(parse_xs_xml_(node.child("xslibs")));
	}

	std::vector<std::size_t> shape = {this->flux_.size()};
	//xt::xarray<double> arr {xt::adapt(this->flux_, shape)};
	xt::xarray<double> arr {0.0, 1.00};
	std::cout << "Xtensor test "<<arr[1] << std::endl;
}


s_xs_ Composition::parse_xs_xml_(pugi::xml_node node) {

	s_xs_ result;
    std::string rxs {node.attribute("typex").value()};
    for (pugi::xml_node tool : node.children("xslib")) {
    	result.xstype = rxs;
    	result.xsname = tool.attribute("name").value();
    	if (rxs == "cs") {
            result.xs_ = get_node_array<double>(tool, "xslib");
    	} else {
    		result.rxs_ = get_node_array<double>(tool, "xslib");
    	}

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
	   this->depcopymap_(this->spectrum_, externcompos.spectrum_);
	   this->depcopymap_(this->flux_, externcompos.flux_);
	   if (!externcompos.namenuclides_.empty() && this->namenuclides_.empty()) {
		   this->namenuclides_.resize(externcompos.namenuclides_.size());
		   std::copy(externcompos.namenuclides_.begin(),externcompos.namenuclides_.end(),
				     this->namenuclides_.begin());
	   }
	   if (!externcompos.namenuclides_.empty() && this->namenuclides_.empty()) {
		   this->conc_.resize(externcompos.conc_.size());
		   std::copy(externcompos.conc_.begin(),externcompos.conc_.end(),
					 this->conc_.begin());
	   }
	   if (!externcompos.xslib_.empty() && this->xslib_.empty()){
		   this->xslib_.resize(externcompos.xslib_.size());
		   std::copy(externcompos.xslib_.begin(),externcompos.xslib_.end(),
		   			 this->xslib_.begin());

	   }
   }
}

void read_reactions_xml() {
	pugi::xml_document doc;
	auto result = doc.load_file(configure::reaction_file.c_str());
	if (!result) {
	    std::cerr << "Error: file not found!" << std::endl;
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
	}


}

}


