#include "nuclide_class.h"
#include "configure.h"
#include "executer.h"
#include "../extern/pugiData/pugixml.h"
#include "materials.h"
std::string XML_INP_PATH = "./Xmls/inpmaterials.xml";
std::string XML_OUT_PATH = "./Xmls/outmaterials.xml";

int main(int arc, char *argv[]) {
    using namespace openbps;
    Timer t;

     int result {parse_command_line(arc, argv)};
     read_input_xml();
     executer::init_solver();

     executer::run_solver();
    //forming from reactions INP
    std::vector<Materials> v = read_materials_from_reactions();
    form_materials_xml(v, XML_INP_PATH);
    //reading v from INP and form to OUT
    v = read_materials_from_inp(XML_INP_PATH);
    form_materials_xml(v, XML_OUT_PATH);
    std :: cout << "Hello world!\n";
    return 0;
}#include "nuclide_class.h"
#include "configure.h"
#include "executer.h"
#include "../extern/pugiData/pugixml.h"
#include "materials.h"
std::string XML_INP_PATH = "./Xmls/inpmaterials.xml";
std::string XML_OUT_PATH = "./Xmls/outmaterials.xml";

int main(int arc, char *argv[]) {
    using namespace openbps;
    Timer t;
//
     int result {parse_command_line(arc, argv)};
     read_input_xml();
     executer::init_solver();

     executer::run_solver();
    //forming from reactions INP
    std::vector<Materials> v = read_materials_from_reactions();
    form_materials_xml(v, XML_INP_PATH);
    //reading v from INP and form to OUT
    v = read_materials_from_inp(XML_INP_PATH);
    form_materials_xml(v, XML_OUT_PATH);
    std :: cout << "Hello world!\n";
    return 0;
}
