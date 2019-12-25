#include "configure.h"
#include "executer.h"
#include "../extern/pugiData/pugixml.h"
#include "chain.h"

int main(int arc, char *argv[]) {
    using namespace openbps;
    Timer t;

     int result {parse_command_line(arc, argv)};
     read_input_xml();
     executer::init_solver();

     executer::run_solver();

    std :: cout << "Hello world!\n";
    return 0;
}
