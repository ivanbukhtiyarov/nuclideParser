#ifndef IBRAE2_NUCLIDE_CLASS_H
#define IBRAE2_NUCLIDE_CLASS_H

#include <iostream>
#include <string.h>
#include <vector>
#include <map>
#include <sstream>
#include "../extern/pugiData/pugixml.h"
#include <chrono>

namespace openbps {

    class Timer {
        using clock_t = std::chrono::high_resolution_clock;
        using microseconds = std::chrono::microseconds;

    public:
        Timer() : start_(clock_t::now()) {}

        ~Timer() {
            const auto finish = clock_t::now();
            const auto us =
                    std::chrono::duration_cast<microseconds>(finish - start_).count();
            std::cout << us << " us" << std::endl;
        }

    private:
        const clock_t::time_point start_;
    };

    struct s_yield {
        double energy;
        std::map<std::string, double> product_data;
    };

    struct s_nfy {
        std::vector<double> energies;
        std::vector<s_yield> yield_arr;
    };

    struct s_decay {
        std::string type;
        std::string target;
        double branching_ratio;
    };

    struct s_reaction {
        std::string type;
        std::string target;
        double q;
    };

    struct s_nuclide {
        std::string name;
        double half_life;
        double decay_energy;
        size_t decay_modes;
        size_t reactions;
        std::vector<s_decay> decay_arr;
        std::vector<s_reaction> reaction_arr;
        s_nfy nfy;
    };

//==============================================================================
// Chain class
//==============================================================================

    class Chain {
    public:
        std::map<std::string, size_t> name_idx;
        std::vector<s_nuclide> nuclides;
        // Constructor
        Chain(){};
        Chain(pugi::xml_node node);
        std::vector<std::pair<int, std::string>> form_idx_name();
        std::vector<std::pair<int, double>> form_idx_lambda();
        std::vector<std::vector<double>> form_idx_decay();
        std::map<std::string, std::vector<std::pair<int, int>>> form_reaction();
        std::map<double, std::vector<std::vector<double>>> form_yield_map();

    private:
        s_decay parse_decay_(pugi::xml_node node);
        s_reaction parse_reaction_(pugi::xml_node node);
        s_yield parse_yield_(pugi::xml_node node);
        s_nfy parse_nfy_(pugi::xml_node node);
        s_nuclide parse_nuclide_(pugi::xml_node node);
    };

    std::string get_node_value(pugi::xml_node node, const char *name);
    bool get_node_value_bool(pugi::xml_node node, const char *name);
    std::vector<double> splitAtof(const std::string &s, char delimiter);
    std::vector<std::string> split(const std::string &s, char delimiter);
    pugi::xml_node read_xml(const std::string& filename);
} // namespace openbps

#endif // IBRAE2_NUCLIDE_CLASS_H
