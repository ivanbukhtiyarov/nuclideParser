#ifndef CHAIN_H
#define CHAIN_H

#include <iostream>
#include <string.h>
#include <vector>
#include <map>
#include <sstream>
#include "../extern/pugiData/pugixml.h"
#include <chrono>
#include "parse.h"
#include "materials.h"
#include "../extern/xtensor/include/xtensor/xarray.hpp"

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
        double d_branching_ratio;
    };

    struct s_reaction {
        std::string type;
        std::string target;
        double q;
    };

    struct s_nuclide {
        std::string name;
        double half_life;
        double d_half_life;
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
        std::pair<std::vector<double>, std::vector<double>>
        get_yield_map_(size_t father, const std::string& daughter);
    private:
        s_decay parse_decay_(pugi::xml_node node);
        s_reaction parse_reaction_(pugi::xml_node node);
        s_yield parse_yield_(pugi::xml_node node);
        s_nfy parse_nfy_(pugi::xml_node node);
        s_nuclide parse_nuclide_(pugi::xml_node node);


    };


    pugi::xml_node read_chain_xml(const std::string& filename);
    xt::xarray<double> form_matrix(Chain& chainer, Materials& mat);
    xt::xarray<double> form_sigp(Chain& chainer, Materials& mat);
    xt::xarray<double> form_dmatrix(Chain& chainer, Materials& mat);
    xt::xarray<double> form_dsigp(Chain& chainer, Materials& mat);
    xt::xarray<double> make_concentration(Chain& chainer, std::vector<std::string>& nameconc,
    		                              std::vector<double>& ro);

} // namespace openbps

#endif // CHAIN_H
