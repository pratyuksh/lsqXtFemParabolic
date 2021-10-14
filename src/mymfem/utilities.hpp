#ifndef MYMFEM_UTILITIES_HPP
#define MYMFEM_UTILITIES_HPP

#include "mfem.hpp"
using namespace mfem;

#include <memory>


//! Deletes memory allocated to SparseMatrix mat
void clear (SparseMatrix* &mat);

//! Defines a function f(x) = 0
void zeroFn(const Vector& x, Vector& f);

//! Returns the upper-triangular, including the diagonal,
//! for the input matrix A
SparseMatrix& getUpperTriangle(const SparseMatrix& A);

//! Returns cell center of element i of given mesh
void getElementCenter(const std::shared_ptr<Mesh>& mesh,
                      int i, Vector& center);

namespace mymfem {

/**
 * @brief Provides functions to locate physical points on a mesh,
 * a faster alternative to the implementation provided by MFEM.
 */
class PointLocator
{
public:
    /**
     * @brief Constructor
     * @param mesh Mesh, passed as pointer
     * @param has_shared_vertices boolean flag when mesh has shared vertices
     */
    PointLocator (Mesh *mesh, bool has_shared_vertices=false)
        : m_hasSharedVertices (has_shared_vertices),
          m_mesh (mesh)
    {
        m_invTr = std::make_unique
                <InverseElementTransformation>();

        if (!m_hasSharedVertices) {
            m_vToEl.reset(m_mesh->GetVertexToElementTable());
        }
        else {
            m_sharedVertices = std::make_unique<Table>();
            getSharedVerticesTable();

            m_vToEl = std::make_unique<Table>();
            getVerticesToElementsTable();
        }
    }

    /**
     * @brief Locates the given physical points on a mesh,
     * speeds up search with an initial guess
     * @param X collection of physical points
     * @param initElId initial guess for index of mesh element
     * @return info flag and collection of reference points
     */
    std::pair <Array<int>, Array<IntegrationPoint>>
    operator() (const DenseMatrix& X, int& initElId) const;

    /**
     * @brief Locates the physical point on a mesh,
     * speeds up search with an initial guess
     * @param x physical point
     * @param initElId initial guess for index of mesh element
     * @return info flag and coordinates of reference point
     */
    std::pair <int, IntegrationPoint> operator()
    (const Vector& x, const int initElId) const;

    /**
     * @brief Generates a table of shared vertices
     */
    void getSharedVerticesTable () const;

    /**
     * @brief Generates vertices to elements connectivity
     */
    void getVerticesToElementsTable () const;

    /**
     * @brief Computes the reference point, given the physical point
     * and mesh element index
     * @param x coordinates of physical point
     * @param elId index of mesh element
     * @return info flag and coordinates of reference point
     */
    inline std::pair <int, IntegrationPoint>
    computeRefPoint (const Vector& x, const int elId) const
    {
        IntegrationPoint ip;
        m_invTr->SetTransformation
                (*m_mesh->GetElementTransformation(elId));
        auto info = m_invTr->Transform(x, ip);
        return {info, ip};
    }

private:
    double m_TOL = 1E-12;
    bool m_hasSharedVertices = false;

    Mesh *m_mesh = nullptr;
    std::unique_ptr<InverseElementTransformation> m_invTr;

    std::unique_ptr<Table> m_vToEl;
    std::unique_ptr<Table> m_sharedVertices;
};

}

#endif // MYMFEM_UTILITIES_HPP
