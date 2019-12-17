#include "reactions.h"
#include <vector>
#include <map>
#include <iostream>
#include <string.h>
#include "configure.h"
#include "nuclide_class.h"
#include "../extern/pugiData/pugixml.h"

namespace openbps {

std::vector<Composition> compositions;
std::map<std::string, int> composmap;
size_t indexall;

Composition::Composition(pugi::xml_node node){

	if (check_for_node(node, "name")) {
	    this->name = node.attribute("name").value();
	}
	if (check_for_node(node, "energy")) {
		for (pugi::xml_node tool : node.children("energy")) {
            size_t cng = std::stoi(tool.attribute("ng").value());
            this->energies_.insert({cng, get_node_array<double>(tool, "energy")})
            this->energy_number++;
		}
	}
	if (check_for_node(node, "namenuclides")) {
	    this->namenuclides_ = get_node_array<std::string>(node, "namenuclides");
	}
	if (check_for_node(node, "flux")) {
	    this->flux_ = get_node_array<std::double>(node, "flux");
	}
	if (check_for_node(node, "spectrum")) {
	    this->spectrum_ = get_node_array<std::double>(node, "spectrum");
	}
	if (check_for_node(node, "conc")) {
		this->conc_ = get_node_array<std::double>(node, "conc");
	}


	pugi::xml_node child_node {node.child("flux")};


}

void Composition::parse_compos_xml_(pugi::xml_node node){

}

s_xs_ Composition::parse_xs_xml_(pugi::xml_node node){

}

void read_reactions_xml(){
	pugi::xml_document doc;
	auto result = doc.load_file(configure::reaction_file.c_str());
	if (!result) {
	    std::cerr << "Error: file not found!" << std::endl;
	}
	pugi::xml_node root_node = doc.child("compositions");

	std::cout << "I' m in reactions.xml parser" << std::endl;

	for (pugi::xml_node tool : root_node.children("composit")) {
		compositions.push_back(tool);
	}


}

}


