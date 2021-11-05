#ifndef MYMFEM_BILINEARFORMS_HPP
#define MYMFEM_BILINEARFORMS_HPP

#include "mfem.hpp"
using namespace mfem;

#include "nested_hierarchy.hpp"
#include "my_bilinearForm_integrators.hpp"


namespace mymfem {

// /**
// * @brief MyBilinearForm class customizes the BilinearForm class provided
// * by MFEM with additional functions.
// */
//class MyBilinearForm : public BilinearForm
//{
//public:
//    //! Constructor with FiniteElementSpace
//    MyBilinearForm(FiniteElementSpace *fes) : BilinearForm(fes) { }
    
//    /**
//     * @brief Adds domain integrators
//     * @param bfi domain bilinear form integrator
//     */
//    void MyAddDomainIntegrator (BilinearFormIntegrator * bfi) {
//        mydbfi.Append (bfi);
//    }

//    /**
//     * @brief Adds interior face integrators
//     * @param bfi interior face bilinear form integrator
//     */
//    void MyAddInteriorFaceIntegrator (BilinearFormIntegrator * bfi) {
//        myifbfi.Append (bfi);
//    }

//    /**
//     * @brief Adds boundary face integrators
//     * @param bfi boundary face bilinear form integrator
//     */
//    void MyAddBoundaryFaceIntegrator (BilinearFormIntegrator * bfi) {
//        mybfbfi.Append (bfi);
//        mybfmarker.Append(nullptr);
//    }

//    /**
//     * @brief Adds boundary face integrators
//     * @param bfi boundary face bilinear form integrator
//     * @param bdr_attr_marker boundary marker
//     */
//    void MyAddBoundaryFaceIntegrator (BilinearFormIntegrator * bfi,
//                                      Array<int> &bdrAttrMarker) {
//        mybfbfi.Append (bfi);
//        mybfmarker.Append(&bdrAttrMarker);
//    }
    
//    /**
//     * @brief Add interior and boundary face integrators
//     * @param bfi face bilinear form integrator
//     */
//    void MyAddFaceIntegrator (BilinearFormIntegrator * bfi) {
//        MyAddInteriorFaceIntegrator(bfi);
//        MyAddBoundaryFaceIntegrator(bfi);
//    }
    
//    /**
//     * @brief Assembles all the integrators
//     * @param skip_zeros Skips zero entries
//     */
//    void MyAssemble(int skip_zeros=1);
    
//    /**
//     * @brief Clears the integrators
//     */
//    void MyClear()
//    {
//        for (int k=0; k < mydbfi.Size(); k++) { delete mydbfi[k]; }
//        for (int k=0; k < myifbfi.Size(); k++) { delete myifbfi[k]; }
//        for (int k=0; k < mybfbfi.Size(); k++) { delete mybfbfi[k]; }
//        for (int k=0; k < myfbfi.Size(); k++) { delete myfbfi[k]; }
//    }

//private:
//    //! domain bilinear form integrator
//    Array<BilinearFormIntegrator *> mydbfi;
//    //! interior face bilinear form integrator
//    Array<BilinearFormIntegrator *> myifbfi;
//    //! boundary face bilinear form integrator
//    Array<BilinearFormIntegrator *> mybfbfi;
//    //! face bilinear form integrator
//    Array<BilinearFormIntegrator *> myfbfi;
//    //! boundary face marker
//    Array<Array<int>*> mybfmarker;
//};


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
    void assemble(int skip_zeros=1);

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
