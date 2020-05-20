#include "../extern/pugiData/pugixml.h"
#include "parse.h"
#include "chain.h"

namespace openbps {
//==============================================================================
// Chain implementation
//==============================================================================

    s_decay Chain::parse_decay_(pugi::xml_node node) {
        s_decay temp;

        temp.type = node.attribute("type").value();
        temp.target = node.attribute("target").value();
        temp.branching_ratio = atof(node.attribute("branching_ratio").value());
        return temp;
    }

    s_reaction Chain::parse_reaction_(pugi::xml_node node) {
        s_reaction temp;

        temp.type = node.attribute("type").value();
        temp.q = atof(node.attribute("Q").value());
        temp.target = node.attribute("target").value();
        return temp;
    }

    s_yield Chain::parse_yield_(pugi::xml_node node) {
        s_yield temp;
        std::vector<std::string> nuclides;
        std::vector<double> numbers;

        temp.energy = atof(node.attribute("energy").value());
        nuclides = split(node.child("products").child_value(), ' ');
        numbers = splitAtof(node.child("data").child_value(), ' ');
        for (std::size_t i = 0; i < nuclides.size(); i++) {
            temp.product_data.insert({nuclides[i], numbers[i]});
        }
        return temp;
    }

    s_nfy Chain::parse_nfy_(pugi::xml_node node) {
        s_nfy temp;

        temp.energies = splitAtof(node.child("energies").child_value(), ' ');
        for (pugi::xml_node tool : node.children("fission_yields")) {
            temp.yield_arr.push_back(parse_yield_(tool));
        }
        return temp;
    }

    s_nuclide Chain::parse_nuclide_(pugi::xml_node node) {
        s_nuclide temp;

        temp.name = node.attribute("name").value();
        temp.decay_modes = atoi(node.attribute("decay_modes").value());
        temp.reactions = atoi(node.attribute("reactions").value());
        temp.half_life = 0;
        temp.half_life = atof(node.attribute("half_life").value());
        std::cout << "Parse start for " << temp.name << std::endl;
        if (temp.decay_modes > 0) {
            temp.decay_energy = atof(node.attribute("decay_energy").value());
            for (pugi::xml_node tool : node.children("decay")) {
                temp.decay_arr.push_back(parse_decay_(tool));
            }
        } else {
            temp.decay_energy = 0;
        }

        if (temp.reactions > 0) {
            for (pugi::xml_node tool : node.children("reaction")) {
                temp.reaction_arr.push_back(parse_reaction_(tool));
            }
        }

        if (node.child("neutron_fission_yields")) {
            temp.nfy = parse_nfy_(node.child("neutron_fission_yields"));
        }
        std::cout << "Parse finish for " << temp.name << std::endl;
        return temp;
    }

    Chain::Chain(pugi::xml_node node) {

        size_t i = 0;
        auto fc = node.first_child();
        std::cout << "First attribute" << fc.first_attribute().value() << "ok\n";
        /*for (pugi::xml_node tool = node.child("nuclide"); tool; tool = tool.next_sibling("nuclide"))
        {

            std::cout << "name " << tool.attribute("name").value();
            std::cout << "reactions " << tool.attribute("reactions").as_int()<< "'\n";
            name_idx.insert({get_node_value(tool, "name"), i});
            nuclides.push_back(parse_nuclide_(tool));
            i++;
        }*/
        for (pugi::xml_node tool : node.children("nuclide")) {
            name_idx.insert({get_node_value(tool, "name"), i});
            nuclides.push_back(parse_nuclide_(tool));
            i++;
        }
        std::cout <<"I'm here!\n ";
    }

} // namespace openbps