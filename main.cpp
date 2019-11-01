#include "nuclide_class.h"

int main(int arc, char* argv[])
{
    using namespace openbps;
    Timer t;
    pugi::xml_document  doc;
    pugi::xml_parse_result result = doc.load_file("./Xmls/chain.xml");
    if (!result) {
        std::cerr << "Error: file not found!" << std::endl;
    }
    
    pugi::xml_node chain_node = doc.child("depletion_chain");

    Chain chain(chain_node);

    std::cout << chain.nuclides["Th232"].reaction_arr[1].q << std::endl;

    return 0;

}
