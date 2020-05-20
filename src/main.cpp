#include "nuclide_class.h"
#include "configure.h"
#include "executer.h"

int main(int arc, char *argv[]) {
    using namespace openbps;
    Timer t;

    int result {parse_command_line(arc, argv)};
    read_input_xml();
    executer::init_solver();

    executer::run_solver();

    return 0;
}
