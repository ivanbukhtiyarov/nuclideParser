#ifndef SRC_MATERIALS_H_
#define SRC_MATERIALS_H_
#include <vector>
#include <map>
#include <iostream>
#include "../extern/pugiData/pugixml.h"

class Materials {
    public:
        Materials();
        void xml_add_material(pugi::xml_node node);
        void check_node(pugi::xml_node node);
        std::string name;
        std::vector<std::string> namenuclides;
        std::vector<double> conc;
        double volume;
        double mass;
        double power;
};

void form_materials_xml(std::vector<Materials> m_arr);
std::vector<Materials> read_materials_from_reactions();
std::vector<Materials> parse_xml_materials();
#endif /* SRC_MATERIALS_H_ */