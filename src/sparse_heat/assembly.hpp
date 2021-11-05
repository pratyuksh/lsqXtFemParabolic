#ifndef SPARSE_HEAT_ASSEMBLY_HPP
#define SPARSE_HEAT_ASSEMBLY_HPP

#include "mfem.hpp"
using namespace mfem;

#include "../mymfem/my_bilinearForms.hpp"


namespace sparseHeat {

/**
 * @brief Mass Integrator; (u, v)
 */
class MassIntegrator
        : public mymfem::BlockBilinearFormIntegrator
{
public:
    //! Default constructor
    explicit MassIntegrator() {}

    //! Assembles the mass integrator on an element
    void assembleElementMatrix
    (const FiniteElement &, ElementTransformation &,
     DenseMatrix &) override;

    //! Assembles the mass integrator between two nested mesh elements
    void assembleElementMatrix
    (const FiniteElement &, ElementTransformation &,
     const FiniteElement &, ElementTransformation &,
     DenseMatrix &) override;
};

}

#endif // SPARSE_HEAT_ASSEMBLY_HPP
