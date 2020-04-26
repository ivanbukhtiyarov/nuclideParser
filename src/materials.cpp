#include "materials.h"
#include "configure.h"
#include "parse.h"

std::string XML_FILE_PATH = "./Xmls/inpmaterials.xml";

Materials::Materials(){
    name = "default material";
    volume =  rand() % 100; ;
    mass =  rand() % 100; ;
    power = 1.0;
}
// void Materials::check_node(pugi::xml_node node) {
// 	if (check_for_node(node, "namenuclides")) {
// 	    this->namenuclides = get_node_array<std::string>(node, "namenuclides");
// 	}

// 	if (check_for_node(node, "conc")) {
// 		this->conc = get_node_array<double>(node, "conc");
// 	}
// }

// std::vector<Materials> parse_xml_materials() {
//     pugi::xml_document doc;
// 	auto result = doc.load_file(configure::reaction_file.c_str());
// 	if (!result) {
// 	    std::cerr << "Error: file not found!" << std::endl;
// 	}
// 	pugi::xml_node root_node = doc.child("compositions");

// 	std::cout << "I' m in material parser" << std::endl;

// 	for (pugi::xml_node tool : root_node.children("composit")) {
// 		compositions.push_back(Composition(tool));
// 		int index = compositions.size() - 1;
// 		if (compositions[index].name == "all") {
// 			indexall = index;
// 		}
// 		composmap.insert({compositions[index].name, index});
// 	}

// 	for (std::vector<Composition>::iterator it = compositions.begin();
// 	     it != compositions.end() ;++it) {
//         it->deploy_all(compositions[indexall]);
//         it->get_reaction();
// 	}
// }

void Materials::xml_add_material(pugi::xml_node node) {
    auto material = node.append_child("material");
    material.append_attribute("name") = name.c_str();
    material.append_attribute("volume") = volume;
    material.append_attribute("mass") = mass;
    material.append_attribute("power") = power;

    auto nameofn = material.append_child("nameofnuclide");
    nameofn.append_child(pugi::node_pcdata).set_value("name nucl");
    auto conc = material.append_child("conc");
    conc.append_child(pugi::node_pcdata).set_value("0 100");
}

void form_materials_xml(std::vector<Materials> m_arr) {
    pugi::xml_document doc;
    auto declarationNode = doc.append_child(pugi::node_declaration);
    declarationNode.append_attribute("version") = "1.0";
    declarationNode.append_attribute("encoding") = "UTF-8";

    auto materials = doc.append_child("materials");
    for(int i = 0; i < m_arr.size(); i++) {
        m_arr[i].xml_add_material(materials);
    }
    bool saveSucceeded = doc.save_file(XML_FILE_PATH.c_str(), PUGIXML_TEXT("  "));
    assert(saveSucceeded);
}   