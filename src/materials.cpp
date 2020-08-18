#include "openbps/materials.h"
#include "openbps/configure.h"
#include "openbps/parse.h"
#include <memory>
#include <algorithm>
#include "../extern/pugiData/pugixml.h"
namespace openbps {

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

    auto nameofn = material.append_child("namenuclides");
    nameofn.append_child(pugi::node_pcdata).set_value(join(namenuclides," ").c_str());
    auto concnod = material.append_child("conc");
    concnod.append_child(pugi::node_pcdata).set_value(joinDouble(conc, " ").c_str());
    std :: cout << "TROUBLE ?" << std ::endl;
    auto dconcnod = material.append_child("dconc");
    std :: cout << "TROUBLE ??" << std ::endl;
    dconcnod.append_child(pugi::node_pcdata).set_value(joinDouble(d_conc, " ").c_str());
    std :: cout << "NO TROUBLE !!" << std ::endl;
}

void Materials::add_nuclide(std::string& extname, double extconc, bool isderiv) {
    auto it = std::find(this->namenuclides.begin(), this->namenuclides.end(), extname);
    if (it == this->namenuclides.end()) {
	this->namenuclides.push_back(extname);
	if (!isderiv) {
	    this->conc.push_back(extconc);
            this->d_conc.push_back(0.0);
	} else {
            std :: cout <<"EVER HERE!"<<std::endl;
	    this->d_conc.push_back(extconc);
            this->conc.push_back(0.0);
	}
    } else {
       auto index = std::distance(this->namenuclides.begin(), it);
       if (!isderiv) {
           this->conc[index] = extconc;
       } else {
	   this->d_conc[index] = extconc;
       }
    }

}

void Materials::bindcomposition(Composition& extcompos){
  this->compos = std::make_shared<Composition>(extcompos);
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

void form_materials_xml(std::vector<Materials> m_arr, std::string xml_path) {
    pugi::xml_document doc;
    auto declarationNode = doc.append_child(pugi::node_declaration);
    declarationNode.append_attribute("version") = "1.0";
    declarationNode.append_attribute("encoding") = "UTF-8";

    auto materials = doc.append_child("materials");
    for(int i = 0; i < m_arr.size(); i++) {
        m_arr[i].xml_add_material(materials);
        std :: cout << "material was added" << std ::endl;
    }
    bool saveSucceeded = doc.save_file(xml_path.c_str(), PUGIXML_TEXT("  "));
    //assert(saveSucceeded);
}

std::vector<Materials> read_materials_from_inp(std::string inp_path) {
    pugi::xml_document doc;
    std::vector<Materials> m_arr;

	auto result = doc.load_file(inp_path.c_str());
	if (!result) {
	    std::cerr << "Error: file not found!" << std::endl;
	}
	pugi::xml_node root_node = doc.child("materials");

	for (pugi::xml_node tool : root_node.children("material")) {
        Materials m;
        m.name = tool.attribute("name").value();
        m.volume = atof(tool.attribute("volume").value());
        m.mass = atof(tool.attribute("mass").value());
        m.power = atof(tool.attribute("power").value());
        auto names = tool.child_value("namenuclides");
        m.namenuclides = split(names,' ');
        auto conc = tool.child_value("conc");
        m.conc = splitAtof(conc, ' ');
        if (check_for_node(tool, "dconc")) {
            auto dconc = tool.child_value("dconc");
            m.d_conc = splitAtof(dconc, ' ');
        }
        m_arr.push_back(m);
	}
    return m_arr;
}

void matchcompositions(std::vector<Composition>& compositions, std::vector<Materials>& materials) {
    for (auto&  mat: materials) {
         for  (auto&  compos:compositions ) {
             if (mat.name==compos.name) {
	              mat.bindcomposition(compos);
                 break;
             }
         }
     }
}

}
