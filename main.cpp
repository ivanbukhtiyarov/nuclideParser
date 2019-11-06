#include "nuclide_class.h"

int main(int arc, char *argv[]) {
    using namespace openbps;
    Timer t;
    pugi::xml_document doc;
    pugi::xml_parse_result result =
            doc.load_file("/Users/ivanbukhtiyarov/ibrae2/Xmls/chain.xml");
    if (!result) {
        std::cerr << "Error: file not found!" << std::endl;
    }

    pugi::xml_node chain_node = doc.child("depletion_chain");

    Chain chain(chain_node);

    //    std::map<std::string, std::vector<std::pair<int, int>>> test;
    //   // chain.form_idx_name();
    //    test = chain.form_reaction();
    //    for (int i = 0; i < test["(n,gamma)"].size(); ++i) {
    //
    //        std::cout <<"first" << test["(n,gamma)"][i].first << "second" <<
    //        test["(n,gamma)"][i].second << std::endl;
    //    }

    return 0;
}
