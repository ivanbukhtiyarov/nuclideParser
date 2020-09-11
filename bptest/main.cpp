#include <iostream>
#include <vector>
#include "nuclide.h"
#include "chain.h"
#include "uncertainty.h"
#include "configure.h"
#include "materials.h"
#include "reactions.h"
#include "functionals.h"
#include "matrix.h"
#include "timeproc.h"
#include <initializer_list>
using namespace std;

const std::string NUCLIDES {"/home/yuri/bptest/nuclides.xml"};
const std::string CHAIN {"/home/yuri/bptest/chain_simple.xml"};

void matrix_test() {
    using namespace openbps;
    std::cout << "MATRIX LOG START\n";
    BaseMatrix<double> m1({{1., 0.}, {1., 1.}});
    BaseMatrix<double> m2({{4., -2.}, {2., 5.}});
    std::cout << m1;
    std::cout << m2;
    m1[1][1] = -3.0;
    double a = m1[1][1];
    m1 = m1 + m2;
    std::cout << m1;
    BaseMatrix<double> m3(2,2);
    BaseMatrix<double> m4 = m1 + m2;
    std::cout << m4;
    m4 = m1 + m3;
    std::cout << m4;
    m4 = m1 | m2;
    std::cout << m4;
    BaseMatrix<double> m5(4,4);
    std::cout << m5;
    m5 = std::move(m2);
    std::cout << m5;
    std::vector<BaseMatrix<double>> vm;
    vm.push_back(BaseMatrix<double>(4,2));
    vm.push_back(std::move(m5));
    vm.push_back(m3);
    vm.push_back(m4);
    for (auto& v : vm[0]) {
        v = 2;
        std::cout << v << " " << std::endl;
    }
    vm.clear();
    std::cout << "MATRIX LOG END\n";
}

int main()
{
    using namespace openbps;
    Timer t;
    read_nuclide_xml(NUCLIDES);
    Chain mychain = read_chain_xml(CHAIN);
    Uncertainty<double> u1(1.0, 1.0);
    Uncertainty<double> u2(2.0, 0.5);
    udouble u3 {u1 + 2.0};
    cout << "Hello World!" << u3.Real() << u3.Dev() << endl;
    cout << "Hello u2: " << u2 << endl;

    std::vector<double> d1 {0.,5.,10.,15.,20.};
    std::vector<double> d2 {0.,10., 20.};
    std::vector<udouble> u4 {{1.,1.},{2.,2.},{3.,3.},{4.,4.}};
    std::vector<std::vector<double>> arr {{1., 0.},{1.,1.0}};
    read_input_xml();
    //matrix_test();
    //std::vector<udouble> u5 {func<udouble>(d1, u4)};
    std::vector<udouble> u5 {collapsing<udouble>(d1, u4, d2)};
    read_materials_from_inp(configure::inmaterials_file);
    form_materials_xml(configure::outmaterials_file);
    read_reactions_xml();
    configure::dumpoutput.resize(10);
    for (auto& d : configure::dumpoutput) {
        d.resize(11);
        for (size_t i = 0; i < 10; i++)
            d[i] = {0.01, 0.0};
    }
    apply_filters(mychain);

    return 0;
}
