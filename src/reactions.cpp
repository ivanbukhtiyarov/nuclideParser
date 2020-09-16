#include "openbps/reactions.h"
#include <vector>
#include <map>
#include <list>
#include <iostream>
#include <memory>
#include <string.h>
#include "openbps/configure.h"
#include "openbps/functionals.h"
#include "openbps/parse.h"
#include "../extern/pugiData/pugixml.h"

namespace openbps {

//==============================================================================
// Global variables
//==============================================================================

std::vector<std::unique_ptr<Composition>> compositions;
std::map<std::string, int> composmap;
size_t indexall;

//==============================================================================
// Composition class implementation
//==============================================================================

Composition::Composition(pugi::xml_node node) {

    if (check_for_node(node, "name")) {
        this->name_ = node.attribute("name").value();
    }
    if (check_for_node(node, "energy")) {
        for (pugi::xml_node tool : node.children("energy")) {
            size_t cng = std::stoi(tool.attribute("ng").value());
            this->energies_.insert({cng, get_node_array<double>(tool,
                                    "energy")});
            this->energy_number_++;
        }
    }
    if (check_for_node(node, "flux")) {
        for (pugi::xml_node tool : node.children("flux")) {
            this->flux_= get_node_array<udouble>(tool, "flux");
        }
    }
    if (check_for_node(node, "dflux")) {
        for (pugi::xml_node tool : node.children("dflux")) {
            std::vector<double> dflux {get_node_array<double>(tool, "dflux")};
            for (size_t j = 0; j < dflux.size() && j < flux_.size(); j++)
                    flux_[j].Adddeviation(dflux[j]);

        }
    }

    if (check_for_node(node, "spectrum")) {
        for (pugi::xml_node tool : node.children("spectrum")) {
            size_t cng = std::stoi(tool.attribute("ng").value());
            this->spectrum_= get_node_array<udouble>(tool, "spectrum");
        }
    }

    if (check_for_node(node, "dspectrum")) {
        for (pugi::xml_node tool : node.children("dspectrum")) {
            std::vector<double> dspectrum {get_node_array<double>(tool,
                                                                  "dspectrum")};
            for (size_t j = 0; j < dspectrum.size() && j< spectrum_.size(); j++)
                    spectrum_[j].Adddeviation(dspectrum[j]);

        }
    }

    if (check_for_node(node, "xslibs")) {
        std::string rxs {node.child("xslibs").attribute("typex").value()};
        for (pugi::xml_node tool : node.child("xslibs").children("xslib")) {
            this->xslib.push_back(parse_xs_xml_(tool, rxs, "xslib"));
        }
        for (pugi::xml_node tool : node.child("xslibs").children("dxslib")) {
            Sxs deriv_rxs {parse_xs_xml_(tool, rxs, "dxslib")};
            for (auto ixs = this->xslib.begin(); ixs != this->xslib.end(); ixs++) {
                if (ixs->xsname == deriv_rxs.xsname &&
                        ixs->xstype == deriv_rxs.xstype) {
                    for (size_t j = 0; j < deriv_rxs.rxs.size() &&
                         j< ixs->rxs.size(); j++)
                            ixs->rxs[j].Adddeviation(deriv_rxs.rxs[j]);
                    for (size_t j = 0; j < deriv_rxs.xs_.size() &&
                         j< ixs->xs_.size(); j++)
                            ixs->xs_[j].Adddeviation(deriv_rxs.xs_[j]);
                }
            }
        }
    }

    this->spectrum_.size() > 0 ?
            this->energy_number_ = this->spectrum_.size() :
            this->energy_number_ = this->flux_.size();
}

//! Auxilary function to copy data from xslib
Sxs parse_xs_xml_
(pugi::xml_node node,const std::string& rxs,const std::string& redex) {
    Sxs result;
    result.xstype = node.attribute("reaction").value();
    result.xsname = node.attribute("name").value();
    if (rxs == "cs") {
        result.xs_ = get_node_array<udouble>(node, redex.c_str());
    } else {
        result.rxs = get_node_array<udouble>(node, redex.c_str());
    }
    return result;
}
//! Auxilary function to copy data from xslib
void Composition::depcopymap_(std::map<size_t, std::vector<double>>& fmap,
                              std::map<size_t, std::vector<double>>& smap) {

    std::map<size_t, std::vector<double>>::iterator it;
    if (!smap.empty()) {
        for ( it= smap.begin(); it != smap.end(); ++it ) {
            if (fmap.find(it->first) == fmap.end()){
                fmap[it->first] = it->second;
            }

        }
    }

}

//! Copy data from composition marked name "all" in xml file
void Composition::deploy_all(Composition& externcompos) {

    if (externcompos.name_ == this->name_) {
        return;
    } else {

        this->depcopymap_(this->energies_, externcompos.energies_);
        if (!externcompos.spectrum_.empty() && this->spectrum_.empty()) {
            this->spectrum_.resize(externcompos.spectrum_.size());
            std::copy(externcompos.spectrum_.begin(),externcompos.spectrum_.end(),
                      this->spectrum_.begin());
        }
        if (!externcompos.flux_.empty() && this->flux_.empty()) {
            this->flux_.resize(externcompos.flux_.size());
            std::copy(externcompos.flux_.begin(),externcompos.flux_.end(),
                      this->flux_.begin());
        }
        if (!externcompos.xslib.empty() && this->xslib.empty()){
            this->xslib.resize(externcompos.xslib.size());
            std::copy(externcompos.xslib.begin(),externcompos.xslib.end(),
                      this->xslib.begin());

        }
    }
}

//! Calculate reaction rate for all reactions in xslib
void  Composition::get_reaction() {
    for (auto ixs = this->xslib.begin(); ixs != this->xslib.end(); ixs++) {
        if (ixs->rxs.empty()) {
            ixs->rxs.resize(1);
            ixs->rxs[0] = 0.0;
            if (!ixs->xs_.empty()) {
                int ng = ixs->xs_.size();
                std::vector<udouble> cuflux(ng, 1.0);
                if (this->flux_.size() > 0) {
                    if (this->flux_.size() == ng) {
                        cuflux = this->flux_;
                    } else {
                        cuflux = collapsing<udouble>(this->energies_[
                                                     this->flux_.size()],
                                                     this->flux_,
                                                     this->energies_[ng]);
                    }
                }
                for (int i = 0; i != ng; i++) {
                    ixs->rxs[0] = ixs->rxs[0] + cuflux[i] * ixs->xs_[i];
                }
            } // if xs_ not empty
        } // if rxs is empty
    }
}
//! Comparator function
bool compfe (std::pair <double, double> a,std::pair <double, double> b) {
    return a.first < b.first;
}

std::pair<std::vector<double>, std::vector<double>>
Composition::get_fluxenergy() {
    std::pair<std::vector<double>, std::vector<double>> result;
    size_t ng;
    double fluxnorm {0.0};
    //std::vector<std::pair<double, double>> forsort_;
    if (spectrum_.empty()) {
        std::for_each(flux_.begin(), flux_.end(), [&] (udouble n) {
            fluxnorm += n.Real();
        });
        for (auto& f: flux_) {
            spectrum_.push_back(f / fluxnorm);
        }
    }
    ng = spectrum_.size();

    result = std::make_pair(energies_[ng], usplit<double>(spectrum_).first);
    return result;
}

//==============================================================================
// Non class methods implementation
//==============================================================================

void read_reactions_xml() {
    using namespace configure;
    pugi::xml_document doc;
    auto result = doc.load_file(reaction_file.c_str());
    if (!result) {
        std::cout << "Warning: file reactions.xml not found!" << std::endl;
        return;
    }
    pugi::xml_node root_node = doc.child("compositions");

    std::cout << "I' m in reactions.xml parser" << std::endl;

    for (pugi::xml_node tool : root_node.children("composit")) {
        compositions.push_back(std::unique_ptr<Composition>(new Composition(tool)));
        int index = compositions.size() - 1;
        if (compositions[index]->Name() == "all") {
            indexall = index;
        }
        composmap.insert({compositions[index]->Name(), index});
    }

    for (size_t i = 0; i < compositions.size(); i++) {
        compositions[i]->deploy_all(*compositions[indexall]);
        compositions[i]->get_reaction();
    }
}

} //namespace openbps
