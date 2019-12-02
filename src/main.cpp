#include "nuclide_class.h"

int main(int arc, char *argv[]) {
    using namespace openbps;
    Timer t;
    std::string input_file{"../Xmls/chain.xml"};
    pugi::xml_node chain_node{read_xml(input_file)};

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

    return 0;
}
