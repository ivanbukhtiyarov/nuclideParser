#include "configure.h"
#include "../extern/pugiData/pugixml.h"
#include "parse.h"
#include <string>

namespace openbps {

//==============================================================================
// Global variable declarations
//==============================================================================

namespace configure {

std::string path_input;            //!< directory where main .xml files resides
std::string path_output {"out.xml"};//!< directory where output files are written
std::string chain_file;            //!< chain-filename.xml
std::string reaction_file;         //!< reaction-filename.xml
int someint {0};                   //!< some int variable
bool somebool {false};             //!< some bool variable
double timestep;                   //!< time of simulation
int numstep;                       //!< number of time step
double epb {1.e-03};               //!< accuracy of calculation
pugi::xml_document docx;           //!< xml document
Mode calcmode {Mode::iteration};    //!< type of solver time depended exponental equation:
//!< 1- direct matrix exponent by default
//!< 2- iteration method by E.F. Seleznev and I.V Chernova 2018
//!< 3- Chebyshev rational approximation by Josey
int order {8};                     //!< CRAM order in {8, 24} by default CRAM16
bool rewrite {true};              //!< whether to rewrite a concentration data by including nuclid from chain
bool outwrite{true};              //!< write calculation result in file
std::vector<std::vector<double>> dumpoutput;     //!< ouput dump
}

//! Parse init line
int parse_command_line(int argc, char* argv[])
{
	for (int i=1; i < argc; ++i) {
        std::string arg {argv[i]};
		if (arg[0] == '-') {
		    if (arg == "-t" || arg == "--test") {
		        configure::someint = 0;
		    	configure::somebool = true;

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
	    configfile = "configure.xml";
	}
	// Parse configure.xml file
	pugi::xml_document doc;
	auto result = doc.load_file(configfile.c_str());
	if (!result) {
	    std::cerr << "Error while processing configure.xml file" ;
	}
	// Get root element
	pugi::xml_node root = doc.document_element();

	  if (check_for_node(root, "chain")) {
	      chain_file = get_node_value(root, "chain");
	  }

	  if (check_for_node(root, "reaction")) {
	      reaction_file = get_node_value(root, "reaction");
	  }

	  if (check_for_node(root, "output")) {
	      path_output = get_node_value(root, "output");
	  }

	  if (check_for_node(root, "numbers")) {
	      numstep = std::stoi(get_node_value(root, "numbers"));
	  }

	  if (check_for_node(root, "timestep")) {
	      timestep = std::stod(get_node_value(root, "timestep"));
	  }

          if (check_for_node(root, "epb")) {
	      epb = std::stod(get_node_value(root, "epb"));
	  }
          if (check_for_node(root, "method")) {
              std::string method = get_node_value(root, "method");
              if (method == "exponent") calcmode = Mode::exponent;
              if (method == "chebyshev") calcmode = Mode::chebyshev;
              if (method == "iteration") calcmode = Mode::iteration;
          }
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
          if (check_for_node(root, "is_outrewrite")) {
              rewrite = get_node_value_bool(root, "is_outrewrite");
          }
}

//! Read input XML files
void read_input_xml()
{
	read_conigure_xml();
}




}





