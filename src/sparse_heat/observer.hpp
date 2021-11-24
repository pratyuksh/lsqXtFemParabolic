#ifndef SPARSE_HEAT_OBSERVER_HPP
#define SPARSE_HEAT_OBSERVER_HPP

#include "mfem.hpp"
using namespace mfem;

#include "../mymfem/base_observer.hpp"


namespace sparseHeat {

/**
 * @brief Observer for the Sparse solver for the Heat equation
 *
 * Provides functions to visualize solutions, compute error
 * and write output to files.
 */
class Observer : public mymfem::BaseObserver
{
public:
    Observer () : BaseObserver () {}

    Observer (const nlohmann::json& config)
              : BaseObserver (config) {}
};

}

#endif // SPARSE_HEAT_OBSERVER_HPP
