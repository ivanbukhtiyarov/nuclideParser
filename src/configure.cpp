#include "configure.h"
#include "../extern/pugiData/pugixml.h"
#include "parse.h"

namespace openbps {

//==============================================================================
// Global variable declarations
//==============================================================================

namespace configure {

std::string path_input;            //!< directory where main .xml files resides
std::string path_output;           //!< directory where output files are written
std::string chain_file;            //!< chain-filename.xml
std::string reaction_file;         //!< reaction-filename.xml
int someint {0};                   //!< some int variable
bool somebool {false};             //!< some bool variable
pugi::xml_document docx;

}

//! Parse init line
int parse_command_line(int argc, char* argv[])
{
	int last_flag = 0;
	for (int i=1; i < argc; ++i) {
        std::string arg {argv[i]};
		if (arg[0] == '-') {
		    if (arg == "-t" || arg == "--test") {
		        configure::someint = 0;
		    	configure::somebool = true;

		      } else {
		        i += 1;
		        configure::someint = 1;
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

	  if (check_for_node(root, "numbers")) {
		  numstep = std::stoi(get_node_value(root, "numbers"));
	  }

	  if (check_for_node(root, "timestep")) {
		  timestep = std::stoi(get_node_value(root, "timestep"));
	  }
}

//! Read input XML files
void read_input_xml()
{
	read_conigure_xml();
}




}




