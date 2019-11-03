//
// Created by Бухтияров  Иван on 02/10/2019.
//

#ifndef IBRAE2_NUCLIDE_CLASS_H
#define IBRAE2_NUCLIDE_CLASS_H

#include <iostream>
#include <string.h>
#include <vector>
#include <map>
#include <sstream>
#include "pugiData/pugixml.h"
#include <chrono>


class Timer
{
    using clock_t = std::chrono::high_resolution_clock;
    using microseconds = std::chrono::microseconds;
public:
    Timer()
            : start_(clock_t::now())
    {
    }

    ~Timer()
    {
        const auto finish = clock_t::now();
        const auto us =
                std::chrono::duration_cast<microseconds>
                        (finish - start_).count();
        std::cout << us << " us" << std::endl;
    }

private:
    const clock_t::time_point start_;
};


typedef struct                      s_yield
{
    double                          energy;
    std::map<std::string,double>    product_data;
}                                   t_yield;

typedef struct                      s_neutron_fission_yields
{
    std::vector<double>             energies;
    std::vector<t_yield>            yield_arr;
}                                   t_nfy;

typedef struct                      s_decay
{
    std::string                     type;
    std::string                     target;
    double                          branching_ratio;
}                                   t_decay;

typedef struct                      s_reaction
{
    std::string                     type;
    std::string                     target;
    double                          q;
}                                   t_reaction;

typedef struct                      s_nuclide
{
    std::string                     name;
    double                          half_life;
    double                          decay_energy;
    size_t                          decay_modes;
    size_t                          reactions;
    std::vector<t_decay>            decay_arr;
    std::vector<t_reaction>         reaction_arr;
    t_nfy                           nfy;
}                                   t_nuclide;

typedef struct                      s_chain
{
    std::map<std::string,t_nuclide> nuclides;
}                                   t_chain;

std::vector<double>                 splitAtof(const std::string& s, char delimiter);
std::vector<std::string>            split(const std::string& s, char delimiter);
t_decay                             parse_decay(pugi::xml_node node);
t_reaction                          parse_reaction(pugi::xml_node node);
t_yield                             parse_yield(pugi::xml_node node);
t_nfy                               parse_nfy(pugi::xml_node node);
t_nuclide                           parse_nuclide(pugi::xml_node node);
t_chain                             parse_chain(pugi::xml_node node);
#endif //IBRAE2_NUCLIDE_CLASS_H
