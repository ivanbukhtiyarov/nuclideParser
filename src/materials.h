#ifndef SRC_MATERIALS_H_
#define SRC_MATERIALS_H_
#include <vector>
#include <map>
#include <iostream>
#include "reactions.h"
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
        Composition* compos;
        void bindcomposition(Composition& extcompos){compos=&extcompos};
        
};

void form_materials_xml(std::vector<Materials> m_arr, std::string xml_path);
std::vector<Materials> read_materials_from_reactions();
std::vector<Materials> read_materials_from_inp(std::string inp_path);
std::vector<Materials> parse_xml_materials();
void matchcompositions(std::vector<Composition>& compositions, std::vector<Materials>& materials);

#endif /* SRC_MATERIALS_H_ */#ifndef SRC_MATERIALS_H_
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

void form_materials_xml(std::vector<Materials> m_arr, std::string xml_path);
std::vector<Materials> read_materials_from_reactions();
std::vector<Materials> read_materials_from_inp(std::string inp_path);
std::vector<Materials> parse_xml_materials();

#endif /* SRC_MATERIALS_H_ */
