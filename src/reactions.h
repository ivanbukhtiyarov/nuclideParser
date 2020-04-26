#ifndef SRC_REACTIONS_H_
#define SRC_REACTIONS_H_
#include <vector>
#include <map>
#include <iostream>
#include "../extern/pugiData/pugixml.h"

namespace openbps {


//?
//struct s_nuclides_ {

//    std::string nuclide_name;
//    double nuclide_conc;
//    int nuclide_ind;
//    double nuclide_mass;

//};

struct s_xs_ {
    std::string xsname;
    std::string xstype;
    std::vector<double> rxs;
    std::vector<double> xs_;
};

class Composition {
public:
    size_t nuclide_number;
    size_t energy_number;
    std::string name;
    std::vector<std::string> namenuclides;
    std::vector<double> conc;
    Composition(){};
    Composition(pugi::xml_node node);
    ~Composition() {
        std::cout << "Composition " + name + " destroyed\n ";};
    void deploy_all(Composition& externcompos);
    void depcopymap_(std::map<size_t, std::vector<double>>& fmap,
    		        std::map<size_t, std::vector<double>>& smap);
    void get_reaction();
    std::pair<std::vector<double>, std::vector<double>> get_fluxenergy();
    std::vector<s_xs_>  xslib;

private:
    //std::vector<s_nuclides_>  nuclides_;

    std::map<size_t, std::vector<double>> energies_;
    std::vector<double> spectrum_;
    std::vector<double> flux_;

    s_xs_ parse_xs_xml_(pugi::xml_node node, std::string rxs);


};



void read_reactions_xml();
extern std::vector<Composition> compositions;
extern std::map<std::string, int> composmap;
extern size_t indexall;
}

#endif /* SRC_REACTIONS_H_ */

