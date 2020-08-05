#ifndef SRC_MATERIALS_H_
#define SRC_MATERIALS_H_
#include <vector>
#include <map>
#include <iostream>
#include <memory>
#include "reactions.h"
#include "../extern/pugiData/pugixml.h"
namespace openbps {
class Materials {
    public:
        Materials();
        void xml_add_material(pugi::xml_node node);
        void check_node(pugi::xml_node node);
        void add_nuclide(std::string& extname, double extconc);
        std::string name;
        std::vector<std::string> namenuclides;
        std::vector<double> conc;
        std::vector<double> d_conc;
        double volume;
        double mass;
        double power;
        std::shared_ptr<Composition> compos;
        void bindcomposition(Composition& extcompos);

};

void form_materials_xml(std::vector<Materials> m_arr, std::string xml_path);
std::vector<Materials> read_materials_from_reactions();
std::vector<Materials> read_materials_from_inp(std::string inp_path);
std::vector<Materials> parse_xml_materials();
void matchcompositions(std::vector<Composition>& compositions, std::vector<Materials>& materials);
}
#endif /* SRC_MATERIALS_H_ */
