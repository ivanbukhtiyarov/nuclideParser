#include "configure.h"
#include <string>
#include <algorithm>
#include <sstream>
#include "pugixml.h"
#include "parse.h"
#include "timeproc.h"
#include "filter.h"

#include "uncertainty.h"
#include "chain.h"
#include "nuclide.h"
#include "materials.h"
#include <cmath>
#include <fstream>
namespace openbps {

//==============================================================================
// Global variable initialization
//==============================================================================

namespace configure {

std::string path_input;           //!< Directory where main .xml files resides
std::string path_output {
    "."};                         //!< Directory where output files are written
std::string chain_file;           //!< Chain-filename.xml
std::string reaction_file;        //!< Reaction-filename.xml
std::string inmaterials_file;     //!< Input materials-filename.xml
std::string outmaterials_file;    //!< Output materials-filename.xml
double timestep;                  //!< Time of simulation
int numstep;                      //!< Number of time step
double epb {1.e-03};              //!< accuracy of calculation
pugi::xml_document docx;          //!< Xml document
Mode calcmode {Mode::iteration};  //!< Type of solver time depended
                                  //!< exponental equation:
                                  //!< 1- direct matrix exponent by default
                                  //!< 2- iteration method by
                                  //!< E.F. Seleznev and I.V Chernova 2018
                                  //!< 3- Chebyshev rational approximation by
                                  //!< Pussa, Josey, etc.
int order {8};                    //!< CRAM order in {8, 24} by default CRAM16
bool rewrite {true};              //!< Whether to rewrite a concentration
                                  //!< data by including nuclid from chain
bool outwrite{true};              //!< Write calculation result in file
std::vector<std::vector<std::array<double, 2>>>
dumpoutput;                       //!< Ouput dump
bool uncertantie_mod;             //!< Calculation mode with uncertanties
                                  //!< taking account
bool decay_extra_out;             //!< Print out more information about
                                  //!< energy decay
std::array<std::string,5>
header_names {"dt",
              "heat",
              "decay-rate",
              "dheat",
              "ddr"};             //!< Name of header in the ouput file

//==============================================================================
// Non class methods implementation
//==============================================================================


//! Parse init line
int parse_command_line(int argc, char* argv[])
{
    for (int i=1; i < argc; ++i) {
        std::string arg {argv[i]};
        if (arg[0] == '-') {
            if (arg == "-t" || arg == "--test") {
                int k {0};
            }
            if (arg == "-o" || arg == "--output") {
                configure::outwrite = true;
            }
        }

    }

    return 0;
}

//! Read configure from XML file
void read_conigure_xml()
{
    using namespace configure;

    std::string configfile;
    if (path_input.length() > 1) {
        configfile = path_input + "configure.xml";
    } else {
        configfile = "/home/yuri/bptest/configure.xml";
    }
    // Parse configure.xml file
    pugi::xml_document doc;
    auto result = doc.load_file(configfile.c_str());
    if (!result) {
        std::cerr << "Error while processing configure.xml file" ;
    }
    // Get root element
    pugi::xml_node root = doc.document_element();
    // Read a name of chain input *.xml file
    if (check_for_node(root, "chain")) {
        chain_file = get_node_value(root, "chain");
    }
    // Read a name of reaction input *.xml file (if presented)
    if (check_for_node(root, "reaction")) {
        reaction_file = get_node_value(root, "reaction");
    }
    // Read a name of materials input *.xml file
    if (check_for_node(root, "inpmaterials")) {
        inmaterials_file = get_node_value(root, "inpmaterials");
    }
    // Read a name of materials output *.xml file
    if (check_for_node(root, "outmaterials")) {
        outmaterials_file = get_node_value(root, "outmaterials");
    }
    // Read an output directory name
    if (check_for_node(root, "output")) {
        path_output = get_node_value(root, "output");
    }
    // Read a number of steps in calculation
    if (check_for_node(root, "numbers")) {
        numstep = std::stoi(get_node_value(root, "numbers"));
    }
    // Read an overal simulation time
    //! TODO: improve time description
    if (check_for_node(root, "timestep")) {
        timestep = std::stod(get_node_value(root, "timestep"));
    }
    // Read a time duration
    //! TODO: improve time description
    if (check_for_node(root, "timerecord")) {
        Timexec duration(root.child("timerecord"));
        timestep = duration.get_seconds();
    }
    // Calculation accuracy (for iteration method)
    if (check_for_node(root, "epb")) {
        epb = std::stod(get_node_value(root, "epb"));
    }
    // Method selection
    if (check_for_node(root, "method")) {
        std::string method = get_node_value(root, "method");
        if (method == "exponent")  calcmode = Mode::exponent;
        if (method == "chebyshev") calcmode = Mode::chebyshev;
        if (method == "iteration") calcmode = Mode::iteration;
    }
    // CRAM order
    if (check_for_node(root, "cram_order")) {
        int ival = std::stoi(get_node_value(root, "cram_order"));
        switch (ival) {
        case 16:
            order = 8;
            break;
        case 48:
            order = 24;
        }
    }
    // Filters
    if (check_for_node(root, "filters")) {
        openbps::read_fitlers_from_xml(root.child("filters"));
    }
    // Simulation keys:
    if (check_for_node(root, "is_outrewrite")) {
        rewrite = get_node_value_bool(root, "is_outrewrite");
    }
    if (check_for_node(root, "decay_print")) {
        outwrite = get_node_value_bool(root, "decay_print");
    }
    if (check_for_node(root, "uncertanties")) {
        uncertantie_mod = get_node_value_bool(root, "uncertanties");
    }
}

} // namespace configure

//! Read input XML files
void read_input_xml()
{
    configure::read_conigure_xml();
}

//! Apply filters to result
void apply_filters(const Chain &chainer) {
    bool isMaterial {false};
    std::array<double, 4> storeval {0., 0., 0., 0.};
    //std::stringstream  output;
    std::ofstream output("/home/yuri/bptest/outlog.csv");
    std::vector<int> xscale;
    std::vector<int> yscale;
    std::vector<int> mainheader;
    std::vector<int> addheader;
    std::vector<std::string> headnames;
    std::vector<std::string> nuclnames;
    xscale.reserve(configure::numstep);
    yscale.reserve(chainer.name_idx.size());
    nuclnames.reserve(chainer.name_idx.size());
    // Applying a time filter to result
    if (timefilter != nullptr) {
        timefilter->apply(configure::timestep / configure::numstep,
                          configure::numstep, xscale);
    } else {
        for (int i = 0; i < configure::numstep; i++)
            xscale.push_back(i);
    }
    // Fill headers and nuclide names array
    std::for_each(chainer.name_idx.begin(),
                  chainer.name_idx.end(),
                  [&nuclnames](std::pair<std::string, size_t> item){
        nuclnames.push_back(item.first);
    });
    std::for_each(configure::header_names.begin(),
                  configure::header_names.end(),
                  [&headnames](std::string item){
        headnames.push_back(item);
    });
    // Search for ordinary string filters
    auto searchexnuclfilter = std::find_if(filters.begin(),
                                           filters.end(),
                                           [] (Filter& f) {
            return f.type == "exnuclide";});
    auto searchnuclfilter = std::find_if(filters.begin(),
                                         filters.end(),
                                         [] (Filter& f) {
            return f.type == "nuclide";});
    auto searchheaderfilter = std::find_if(filters.begin(),
                                           filters.end(),
                                           [] (Filter& f) {
            return f.type == "header";});
    // If nuclide filter is present -> apply it
    if (searchexnuclfilter != filters.end()) {
        searchexnuclfilter->apply(nuclnames, yscale);
    }
    // If header filter is present -> apply it
    if (searchheaderfilter != filters.end()) {
        mainheader.push_back(0);
        searchheaderfilter->apply(headnames, mainheader);
    } else {
        for (int i = 0; i < configure::header_names.size(); i++)
            mainheader.push_back(i);
    }
    // Addition to header if nuclide concentration specifier is present
    if (searchnuclfilter != filters.end()) {
        searchnuclfilter->apply(nuclnames, addheader);
    }
    for (size_t k = 0; k < mainheader.size(); k++) {
        switch (mainheader[k]) {
        case 0:
            output << "dt" << ";";
            break;
        case 1:
            output << "Act, sec-1" << ";";
            break;
        case 2:
            output << "Q, Mev" << ";";
            break;
        case 3:
            if (configure::uncertantie_mod)
                output << "dAct, sec-1" << ";";
            break;
        case 4:
            if (configure::uncertantie_mod)
                output << "dQ, Mev" << ";";
            break;
        }
    }
    for (auto& k : addheader)
        output << nuclnames[k] << ";";
    output << std::endl;
    for (auto& m : materials) {
        // If materials fitler is present
        // then applying it
        if (materialfilter != nullptr) {
            materialfilter->apply(m->Name(), isMaterial);
        } else {
            isMaterial = true;
        }
        if (isMaterial) {
            output << m->Name() << ";" << std::endl;
            for (auto t : xscale) {
                output << configure :: timestep / configure :: numstep *
                          (t + 1) << ";";
                for (size_t i = 0; i < nuclnames.size(); i++) {
                    if (std::find(yscale.begin(), yscale.end(), i) ==
                            yscale.end()) {
                        int ijk {chainer.name_idx.at(
                                        nuclnames[i])};
                        udouble uhalflife {nuclides[chainer.name_idx.at(
                                        nuclnames[i])]->half_life};
                        for (size_t k = 0; k < mainheader.size(); k++) {
                            if (mainheader[k] == 1 &&
                                uhalflife.Real() > 0)//Decay-rate
                                storeval[0] +=
                                        log(2.0) /
                                        uhalflife.Real()
                                        * configure::dumpoutput[t][i][0] *
                                        1.e+24;
                            if (mainheader[k] == 2 &&
                                    uhalflife.Real() > 0) //Decay heat
                                storeval[1] +=
                                        log(2.0) /
                                        uhalflife.Real() *
                                        nuclides[ijk]->decay_energy
                                        * configure::dumpoutput[t][i][0] *
                                        1.e+24 * 1.e-6;
                            if (mainheader[k] == 3 &&
                                    uhalflife.Real() > 0) //d(Decay-rate)
                                storeval[2] +=
                                        (log(2.0) /
                                        uhalflife).Dev() *
                                        configure::dumpoutput[t][i][0] *
                                        1.e+24 +
                                        log(2.0) /
                                        uhalflife.Real()
                                        * configure::dumpoutput[t][i][1] *
                                        1.e+24;
                            if (mainheader[k] == 4 &&
                                    uhalflife.Real() > 0) //d(Decay heat)
                                storeval[3] +=
                                        ((log(2.0) /
                                          uhalflife
                                        ).Dev() *
                                        configure::dumpoutput[t][i][0] +
                                        log(2.0) /
                                        uhalflife.Real()
                                        * configure::dumpoutput[t][i][1]) *
                                        nuclides[ijk]->decay_energy *
                                        1.e+24 * 1.e-6;
                        } // for with if
                    }
                } // for nuclides
                for (size_t k = 1; k < mainheader.size(); k++)
                    output << storeval[mainheader[k] - 1] << ";";

                for (auto& k : addheader)
                output << configure::dumpoutput[t][k][0] << ";";
                output << std::endl;
            } // for time
        } // if
        for (size_t k = 1; k < mainheader.size(); k++)
            storeval[mainheader[k] - 1] = 0.0;
    } // materials
    output.close();
    /*std::ofstream myfile("/home/yuri/bptest/outlog1.csv");
    if (myfile.is_open()) {
        myfile << output.str();
        myfile.close();
    }*/
}
} // namespace openbps





