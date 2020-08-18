#ifndef SRC_CONFIGURE_H_
#define SRC_CONFIGURE_H_

#include "../extern/pugiData/pugixml.h"
#include <sstream>
#include <vector>
namespace openbps {


//==============================================================================
// Global variable declarations
//==============================================================================

namespace configure {

enum class Mode
  {exponent, //simple exponent matrix
   iteration,//iterative alogrithm
   chebyshev // Chebyshev rational approximation method
  };
extern std::string path_input;            //!< directory where main .xml files resides
extern std::string path_output;           //!< directory where output files are written
extern std::string chain_file;            //!< chain-filename.xml
extern std::string reaction_file;         //!< reaction-filename.xml
extern std::string inmaterials_file;      //!< input materials-filename.xml
extern std::string outmaterials_file;     //!< output materials-filename.xml
extern int someint;                       //!< some int variable
extern int numstep;                       //!< number of substep per one time step
extern double timestep ;                  //!< length of time interval for decesion
extern double epb;                        //!< accuracy of calculation
extern bool somebool;                     //!< some bool variable
extern pugi::xml_document docx;           //!< file handle for chain.xml
extern Mode calcmode;                     //!< mode of calculation
extern int order;                         //!< CRAM order in {8, 24}
extern bool rewrite;                      //!< whether to rewrite a concentration data by including nuclid from chain
extern bool outwrite;                     //!< write calculation result in file
extern std::vector<std::vector<double>> 
       dumpoutput;                        //!< ouput dump
extern bool uncertantie_mod;              //!< calculation mode with uncertanties taking account
extern bool decay_extra_out;              //!< print out more information about energy decay
}

//! Parse init line
int parse_command_line(int argc, char* argv[]);

//! Read configure from XML file
void read_conigure_xml();

//! Read input XML files
void read_input_xml();

}

#endif /* SRC_CONFIGURE_H_ */
