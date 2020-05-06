#include "materials.h"
#include "configure.h"
#include "parse.h"

std::string XML_FILE_PATH = "./Xmls/inpmaterials.xml";
using namespace openbps;

Materials::Materials(){
    name = "default material";
    volume =  rand() % 100; ;
    mass =  rand() % 100; ;
    power = 1.0;
} 


void Materials::xml_add_material(pugi::xml_node node) {
    auto material = node.append_child("material");
    material.append_attribute("name") = name.c_str();
    material.append_attribute("volume") = volume;
    material.append_attribute("mass") = mass;
    material.append_attribute("power") = power;

    auto nameofn = material.append_child("nameofnuclide");
    nameofn.append_child(pugi::node_pcdata).set_value(join(namenuclides," ").c_str());
    auto concnod = material.append_child("conc");
    concnod.append_child(pugi::node_pcdata).set_value(joinDouble(conc, " ").c_str());
}

std::vector<Materials> read_materials_from_reactions() {
    pugi::xml_document doc;
    std::vector<Materials> m_arr;

	auto result = doc.load_file(configure::reaction_file.c_str());
	if (!result) {
	    std::cerr << "Error: file not found!" << std::endl;
	}
	pugi::xml_node root_node = doc.child("compositions");

	std::cout << "I' m in reactions.xml parser (FOR MATERIALS)" << std::endl;

	for (pugi::xml_node tool : root_node.children("composit")) {
        Materials m;
        std::cout << tool.first_child().name() << std::endl;
        auto names = tool.child_value("namenuclides");
        m.namenuclides = split(names,' ');
        auto conc = tool.child_value("conc");
        m.conc = splitAtof(conc, ' ');
        m_arr.push_back(m);
	}
    return m_arr;
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