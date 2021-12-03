#ifndef MYMFEM_BILINEARFORMS_HPP
#define MYMFEM_BILINEARFORMS_HPP

#include "mfem.hpp"

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
    mfem::Array<int> getBlockOffsets() {
        return m_blockOffsets;
    }

    //! Returns the block matrix
    std::shared_ptr<mfem::BlockMatrix> getBlockMatrix() {
        return m_blockMatrix;
    }

private:
    const std::shared_ptr<NestedFEHierarchy> m_nestedFEHierarchy;

    int m_numLevels;
    mfem::Array<int> m_blockOffsets;

    //! Block Matrix, sparse
    std::shared_ptr<mfem::BlockMatrix> m_blockMatrix;

    //! domain block bilinear form integrators
    std::vector<std::shared_ptr<BlockBilinearFormIntegrator>> mydbfi;
};


/**
 * @brief BlockMixedBilinearForm class provides functions to
 * assemble mixed bilinear forms between different mesh levels
 * in a nested hierarchy of finite element spaces
 */
class BlockMixedBilinearForm
{
public:
    //! Constructor with test and trial
    //! nested finite element hierarchies passed as arguments
    BlockMixedBilinearForm
    (const std::shared_ptr<NestedFEHierarchy>& trialNestedFEHierarchy,
     const std::shared_ptr<NestedFEHierarchy>& testNestedFEHierarchy);

    /**
     * @brief Adds domain integrators
     * @param bfi domain bilinear form integrator
     */
    void addDomainIntegrator
    (std::shared_ptr<BlockMixedBilinearFormIntegrator>& bfi) {
        mydbfi.push_back(bfi);
    }

    /**
     * @brief Assembles all the integrators
     * @param skip_zeros Skips zero entries
     */
    void assemble(int skip_zeros=1);

    /**
     * @brief Allocates memory for the block matrix member variable
     */
    void allocateBlockMatrix();

    //! Returns the row block offsets
    mfem::Array<int> getRowBlockOffsets() {
        return m_rowBlockOffsets;
    }

    //! Returns the column block offsets
    mfem::Array<int> getColBlockOffsets() {
        return m_colBlockOffsets;
    }

    //! Returns the block matrix
    std::shared_ptr<mfem::BlockMatrix> getBlockMatrix() {
        return m_blockMatrix;
    }

private:
    const std::shared_ptr<NestedFEHierarchy> m_trialNestedFEHierarchy;
    const std::shared_ptr<NestedFEHierarchy> m_testNestedFEHierarchy;

    int m_numLevels;
    mfem::Array<int> m_rowBlockOffsets, m_colBlockOffsets;

    //! Block Matrix, sparse
    std::shared_ptr<mfem::BlockMatrix> m_blockMatrix;

    //! domain block bilinear form integrators
    std::vector<std::shared_ptr<BlockMixedBilinearFormIntegrator>> mydbfi;
};

}

#endif // MYMFEM_BILINEARFORMS_HPP
