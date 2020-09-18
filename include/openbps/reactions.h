#ifndef SRC_REACTIONS_H_
#define SRC_REACTIONS_H_
#include <vector>
#include <list>
#include <map>
#include <iostream>
#include <algorithm>
#include <memory>
#include "../extern/pugiData/pugixml.h"
#include "uncertainty.h"

namespace openbps {

//==============================================================================
// Global variables
//==============================================================================

class BaseCompostion;
class Composition;

extern std::vector<
std::unique_ptr<Composition>> compositions;   //!< vector with all compositions
extern int indexall;                          //!< index of "all" record with
                                              //!< common shared data
extern std::map<std::string, int> composmap;  //!< composition map with name:int
                                              //!< values

//==============================================================================
// Struct cross-section description
//==============================================================================

struct Sxs {
    std::string xsname;      //!< Name of cross-section/reaction
    std::string xstype;      //!< Cross-section/reaction type
    std::vector<udouble> rxs;//!< Reactions vector by energy group
    std::vector<udouble> xs_;//!< Cross-sections values by energy group

};

//==============================================================================
// Compositions classes descriptions
//==============================================================================
class BasicComposition {
public:
    //--------------------------------------------------------------------------
    // Constructors, destructors, factory functions
    BasicComposition() {std::cout << "basic created"<<std::endl;}
    BasicComposition(std::string name, size_t nuclidnumber, size_t energynumber):
                     name_{name}, nuclide_number_{nuclidnumber},
                     energy_number_{energynumber} {}
    virtual ~BasicComposition() = default;
    //--------------------------------------------------------------------------
    // Methods
    //! Get name
    //!
    //! \return Composition name
    std::string Name() {return name_;}
    //! Get nuclid number
    //!
    //! \return number of nuclides
    size_t NuclidNumber() {return nuclide_number_;}
    //! Get a number of energy discretezation interval
    //!
    //! \return number of energy points
    size_t EnergyNumber() {return energy_number_;}

protected:
    //--------------------------------------------------------------------------
    // Attributes
    size_t nuclide_number_; //!< Nuclide number
    size_t energy_number_;  //!< Energy points number
    std::string name_;      //!< Composition Name
};

class Composition : public BasicComposition {
public:
    //--------------------------------------------------------------------------
    // Constructors, destructors, factory functions
    Composition(Composition& tmp) = delete;
    Composition& operator =(const Composition& tmp) = delete;
    Composition(Composition&& tmp) = default;
    Composition&& operator =(Composition&& tmp) = delete;
    Composition(std::string name, size_t nuclidnumber, size_t energynumber):
                     BasicComposition(name, nuclidnumber, energynumber) {}

    Composition(pugi::xml_node node);
    //--------------------------------------------------------------------------
    // Methods
    //! Copy data from composition marked name "all" in xml file
    //!
    //! \param[in] externcompos external composition
    void deploy_all(Composition &externcompos);
    //! Calculate reaction rate for all reactions in xslib
    void get_reaction();
    //! Get spectrum energy distribution
    //!
    //! \return pair with energies and spectrum distr
    std::pair<std::vector<double>, std::vector<double>> get_fluxenergy();
    //--------------------------------------------------------------------------
    // Attributes
    std::vector<Sxs> xslib; //!< Cross-section/reactions data

private:
    //--------------------------------------------------------------------------
    // Methods
    //! Auxilary function to copy data from xslib
    //!
    //! \param[in] fmap external composition data
    //! \param[out] smap the copied data
    void depcopymap_(std::map<size_t, std::vector<double>>& fmap,
                     std::map<size_t, std::vector<double>>& smap);

    //--------------------------------------------------------------------------
    // Attributes
    std::map<size_t, std::vector<double>> energies_; //!< Energies
                                                     //!< discretization
    std::vector<udouble> spectrum_;                  //!< Energy spectrum
    std::vector<udouble> flux_;                      //!< Energy flux

};

//==============================================================================
// Non class methods
//==============================================================================
//! Auxilary function to copy data from xslib
//!
//! \param[in] node xml node element with xslib name
//! \param[in] rxs the sign of reactions (rxs)/ cross-section(xs) data
Sxs parse_xs_xml_
(pugi::xml_node node, const std::string& rxs,const std::string& redex);

//! Read compositions from xml file
void read_reactions_xml();

}

#endif /* SRC_REACTIONS_H_ */

