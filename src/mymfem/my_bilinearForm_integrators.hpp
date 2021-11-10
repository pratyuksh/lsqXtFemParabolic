#ifndef MYMFEM_BILINEARFORM_INTEGRATORS_HPP
#define MYMFEM_BILINEARFORM_INTEGRATORS_HPP

#include "mfem.hpp"
using namespace mfem;

#include <fmt/format.h>


namespace mymfem {

/**
 * @brief BlockBilinearForm class provides functions to
 * assemble bilinear forms between different mesh levels
 * in a nested hierarchy of finite element spaces
 */
class BlockBilinearFormIntegrator
{
public:
    //! Default destructor
    virtual ~ BlockBilinearFormIntegrator() = default;

    //! Assembles the integrator over a spatial element
    virtual void assembleElementMatrix
    (const FiniteElement &fe,
     ElementTransformation& elTrans,
     DenseMatrix& elmat) = 0;

    //! Assembles the integrator
    //! over two nested coarse and fine mesh elements.
    //! Coarse element is for the trial functions
    //! and fine element is for the test function
    virtual void assembleElementMatrix
    (const FiniteElement &testFeFine,
     ElementTransformation& testElTransFine,
     const FiniteElement &trialFeCoarse,
     ElementTransformation& trialElTransCoarse,
     DenseMatrix& elmat) = 0;
};

/**
 * @brief BlockMixedBilinearForm class provides functions to
 * assemble mixed bilinear forms between different mesh levels
 * in a nested hierarchy of finite element spaces
 */
class BlockMixedBilinearFormIntegrator
{
public:
    //! Default destructor
    virtual ~ BlockMixedBilinearFormIntegrator() = default;

    //! Assembles the integrator over a spatial element
    virtual void assembleElementMatrix
    (const FiniteElement &testFe,
     const FiniteElement &trialFe,
     ElementTransformation& elTrans,
     DenseMatrix& elmat) = 0;

    //! Assembles the integrator
    //! over two nested coarse and fine mesh elements.
    //! Coarse element is for the trial functions
    //! and fine element is for the test function
    virtual void assembleElementMatrix
    (const FiniteElement &testFeFine,
     ElementTransformation& testElTransFine,
     const FiniteElement &trialFeCoarse,
     ElementTransformation& trialElTransCoarse,
     DenseMatrix& elmat) = 0;

    //! Assembles the integrator
    //! over two nested coarse and fine mesh elements.
    //! Fine element is for the trial functions
    //! and coarse element is for the test function
    virtual void assembleElementMatrix2
    (const FiniteElement &testFeCoarse,
     ElementTransformation& testElTransCoarse,
     const FiniteElement &trialFeFine,
     ElementTransformation& trialElTransFine,
     DenseMatrix& elmat) = 0;
};

}

#endif // MYMFEM_BILINEARFORM_INTEGRATORS_HPP
