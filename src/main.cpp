#include "nuclide_class.h"
#include "configure.h"
#include "executer.h"
#include "../extern/pugiData/pugixml.h"
#include "materials.h"
 //std::string XML_OUT_PATH = "./Xmls/outmaterials.xml";

int main(int arc, char * argv[]) {
    using namespace openbps;
    Timer t;

    int result {
        parse_command_line(arc, argv)
    };
    read_input_xml();
    executer::init_solver();

    executer::run_solver();
    //to create new inp_materials
    // std::vector<Materials> v = read_materials_from_reactions();
    // form_materials_xml(v, configure::inp_materials_file);
    //forming from reactions INP
    auto w = read_materials_from_inp(configure::inp_materials_file);
    //form_materials_xml(v, XML_OUT_PATH);
    std::cout << "Hello world!\n";
    return 0;
}