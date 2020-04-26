#include "nuclide_class.h"
#include "configure.h"
#include "executer.h"
#include "../extern/pugiData/pugixml.h"
#include "materials.h"

int main(int arc, char *argv[]) {
    using namespace openbps;
    Timer t;

     int result {parse_command_line(arc, argv)};
     read_input_xml();
     executer::init_solver();

     executer::run_solver();
//some random materials
    std::vector<Materials> v;
    for(int i = 0; i<= 200; i++) {
        Materials m;
        v.push_back(m);
    }
    form_materials_xml(v);
//some random materials
    std :: cout << "Hello world!\n";
    return 0;
}
