#ifndef SRC_CONFIGURE_H_
#define SRC_CONFIGURE_H_

#include "../extern/pugiData/pugixml.h"

namespace openbps {

//==============================================================================
// Global variable declarations
//==============================================================================

namespace configure {

extern std::string path_input;            //!< directory where main .xml files resides
extern std::string path_output;           //!< directory where output files are written
extern std::string chain_file;            //!< chain-filename.xml
extern std::string reaction_file;         //!< reaction-filename.xml
extern int someint;                       //!< some int variable
extern int numstep;                       //!< number of substep per one time step
extern double timestep ;                  //!< length of time interval for decesion
extern bool somebool;                     //!< some bool variable
extern pugi::xml_document docx;           //!< file handle for chain.xml


}

//! Parse init line
int parse_command_line(int argc, char* argv[]);

//! Read configure from XML file
void read_conigure_xml();

//! Read input XML files
void read_input_xml();

}

#endif /* SRC_CONFIGURE_H_ */
