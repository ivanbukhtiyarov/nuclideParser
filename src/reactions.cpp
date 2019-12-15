#include "reactions.h"
#include <vector>
#include <map>
#include <iostream>
#include <string.h>

namespace openbps {

Composition::Composition(pugi::xml_node node){

}

s_nuclides_ Composition::parse_nuclides_xml_(pugi::xml_node node){

}
    s_xs_       Composition::parse_xs_xml_(pugi::xml_node node){

}

void read_reactions_xml(){
	//using namespace openbps;

	std::cout << "I' m in reactions.xml parser" << std::endl;

}

}


