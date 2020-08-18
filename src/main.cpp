#include "openbps/nuclide_class.h"
#include "openbps/configure.h"
#include "openbps/executer.h"

int main(int arc, char *argv[]) {
    using namespace openbps;
    Timer t;

    int result {parse_command_line(arc, argv)};
    read_input_xml();
    executer::init_solver();

    executer::run_solver();

    return 0;
}
