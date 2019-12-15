#include "executer.h"
#include "configure.h"
#include "nuclide_class.h"
#include "functionals.h"

namespace openbps {

namespace executer {

void run_solver() {

	std::vector<double> xval = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
	std::vector<double> yval = {10.0, 22.0, 32.0, 33.0, 14.0, 85.0};
	std::vector<double> xtarget = {0.0, 6.0, 8.2};

//	std::vector<double> ytarget {collapsing(xval, yval, xtarget)};
	std::cout << "We star execution\n";

}

void init_solver() {

	//std::string input_file{"./Xmls/chain.xml"};
	pugi::xml_node chain_node{read_chain_xml(configure::chain_file)};

	Chain chain(chain_node);
	std::cout << "Simple reaction\n";
	    auto react_map = chain.form_reaction();
	    for(int i = 0; i < react_map["(n,gamma)"].size(); i++)
	    {
	        std::cout << react_map["(n,gamma)"][i].first << "  " << react_map["(n,gamma)"][i].second << "  "<< 1 << "  " <<std::endl;
	    }

	    std::cout << "NFY reaction\n";
	    auto map = chain.form_yield_map();
	    std::cout << "Energies of reaction\n";
	    for(auto iter = map.begin() ; iter != map.end() ; iter++)
	    {
	        std::cout << "Energy is " << iter->first << std::endl;
	    }
	    for(int i = 0; i < map[0.0253].size(); i++)
	    {
	        std::cout << map[0.0253][i][0] << "  " << map[0.0253][i][1] << "  "<< map[0.0253][i][2] << "  " <<std::endl;
	    }

}


}

}


