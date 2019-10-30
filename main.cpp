
#include "nuclide_class.h"

int main()
{
    Timer t;
    pugi::xml_document  doc;
    t_decay             decay;
    t_yield             yield;
    t_nfy               nfy;
    t_nuclide           nuclide;
    pugi::xml_parse_result result = doc.load_file("/Users/ivanbukhtiyarov/Desktop/c++/ibrae2/Xmls/chain.xml");
    t_chain chain;
    pugi::xml_node chain_node = doc.child("depletion_chain");
//    pugi::xml_node decay_node = doc.child("depletion_chain").child("nuclide").child("decay");
//    pugi::xml_node yield_node = doc.child("depletion_chain").child("nuclide").child("neutron_fission_yields").child("fission_yields");
//    pugi::xml_node nfy_node = doc.child("depletion_chain").child("nuclide").child("neutron_fission_yields");
//    pugi::xml_node nuclide_node = doc.child("depletion_chain").child("nuclide");

    chain = parse_chain(chain_node);
//    decay = parse_decay(decay_node);
//    yield = parse_yield(yield_node);
//    nfy = parse_nfy(nfy_node);
//    nuclide = parse_nuclide(nuclide_node);
    std::cout << chain.nuclides["Rg272"].decay_energy << std::endl;

    return 0;

    // tag::basic[]
    // end::contents[]
}