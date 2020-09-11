#include "matrix.h"
#include <vector>
#include <memory>
#include <initializer_list>
#include <cmath>
#include "uncertainty.h"
#include "nuclide.h"
#include "chain.h"
#include "reactions.h"

namespace openbps {

//==============================================================================
// Decay matrix implementation
//==============================================================================

//! Form a decay nuclide matrix
void DecayMatrix::form_matrixreal(Chain& chain) {
    size_t NN {chain.name_idx.size()}; //!< Nuclide number
    size_t i {0};
    size_t k;
    udouble decay_ {0.0, 0.0};
    //Run over all nuclides
    for (auto it = chain.name_idx.begin(); it != chain.name_idx.end(); it++) {
        i = it->second; // from nuclide
        if (nuclides[i]->half_life.Real() > 0) {
            decay_ = (log(2.0) / nuclides[i]->half_life);
            this->data_[i][i] = -decay_.Real();
            // Cycle over all decay modes
            for (auto& v : nuclides[i]->get_decaybr()) {
                if (v.first != std::string("Nothing")) {
                    // To nuclide
                    k = chain.name_idx[v.first];
                    this->data_[k][i] +=
                            decay_.Real() * v.second.Real();
                }
            }

        }
    }

}

//! Form a deviation decay nuclide matrix
void DecayMatrix::form_matrixdev(Chain& chain) {
    size_t NN {chain.name_idx.size()}; //!< Nuclide number
    size_t i {0};
    size_t k;
    udouble decay_ {0.0, 0.0};
    //Run over all nuclides
    for (auto it = chain.name_idx.begin(); it != chain.name_idx.end(); it++) {
        i = it->second; // from nuclide
        if (nuclides[i]->half_life.Real() > 0) {
            decay_ = (log(2.0) / nuclides[i]->half_life);
            this->data_[i][i] = -decay_.Dev();
            // Cycle over all decay modes
            for (auto& v : nuclides[i]->get_decaybr()) {
                if (v.first != std::string("Nothing")) {
                    // To nuclide
                    k = chain.name_idx[v.first];
                    this->data_[k][i] +=
                            (decay_ * v.second).Dev();
                }
            }

        }
    }

}

} //namespace openbps
