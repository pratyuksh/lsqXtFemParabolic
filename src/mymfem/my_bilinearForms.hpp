#ifndef MYMFEM_BILINEARFORMS_HPP
#define MYMFEM_BILINEARFORMS_HPP

#include "mfem.hpp"
using namespace mfem;

#include "nested_hierarchy.hpp"
#include "my_bilinearForm_integrators.hpp"


namespace mymfem {

/**
 * @brief BlockBilinearForm class provides functions to
 * assemble bilinear forms between different mesh levels
 * in a nested hierarchy of finite element spaces
 */
class BlockBilinearForm
{
public:
    //! Constructor with nested finite element hierarchy as argument
    BlockBilinearForm(const std::shared_ptr<NestedFEHierarchy>&
                      nestedFEHierarchy);

    /**
     * @brief Adds domain integrators
     * @param bfi domain bilinear form integrator
     */
    void addDomainIntegrator
    (std::shared_ptr<BlockBilinearFormIntegrator>& bfi) {
        mydbfi.push_back(bfi);
    }

    /**
     * @brief Assembles all the integrators
     * @param skip_zeros Skips zero entries
     */
    void assemble(bool symmetric=true, int skip_zeros=1);

    /**
     * @brief Allocates memory for the block matrix member variable
     */
    void allocateBlockMatrix();

    //! Returns the block offsets
    Array<int> getBlockOffsets() {
        return m_blockOffsets;
    }

    //! Returns the block matrix
    std::shared_ptr<BlockMatrix> getBlockMatrix() {
        return m_blockMatrix;
    }

private:
    const std::shared_ptr<NestedFEHierarchy> m_nestedFEHierarchy;

    int m_numLevels;
    Array<int> m_blockOffsets;

    //! Block Matrix, sparse
    std::shared_ptr<BlockMatrix> m_blockMatrix;

    //! domain block bilinear form integrators
    std::vector<std::shared_ptr<BlockBilinearFormIntegrator>> mydbfi;
};


}

#endif // MYMFEM_BILINEARFORMS_HPP
