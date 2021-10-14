#ifndef MYMFEM_MESH_HIERARCHY_HPP
#define MYMFEM_MESH_HIERARCHY_HPP

#include "mfem.hpp"
using namespace mfem;

#include <vector>
#include <set>
#include <memory>

#include "utilities.hpp"


typedef struct MeshHierarchyTable
{
public:
    //! Constructor with number of coarse mesh elements
    MeshHierarchyTable(int numElsCoarse) : m_numRows(numElsCoarse) {
        m_data.resize(m_numRows);
    }

    //! Adds an element in the fine mesh to the hierarchy of an element in
    //! the coars mesh
    void add (int elIdCoarse, int elIdFine) {
        m_data[elIdCoarse].insert(elIdFine);
    }
   
    //! Returns the child indices for the coarse mesh element with
    //! index elIdCoarse
    std::set<int>& operator() (int elIdCoarse) {
        return m_data[elIdCoarse];
    }

    //! Returns the number of rows
    int getNumRows() const {
        return m_numRows;
    }

    //! Returns the number of children for the coarse mesh element with
    //! index elIdCoarse
    int getNumChildren(int elIdCoarse) const {
        return m_data[elIdCoarse].size();
    }

   void print() {
       for (int i=0; i<m_numRows; i++) {
           for (auto it=m_data[i].begin(); it!= m_data[i].end(); it++) {
               std::cout << *it << "\t";
           }
           std::cout << "\n";
       }
   }

private:
   int m_numRows; // number of coarse mesh elements
   std::vector<std::set<int>> m_data; // stores the mesh hierarchy
} MeshHierarchyTable;


namespace mymfem
{
/**
 * @brief Provides routines to build the hierarchy
 * between a given coarse and fine mesh
 */
class MeshHierarchy {
public:
    //! Constructor with coarse and fine mesh passed as arguments
    MeshHierarchy(const std::shared_ptr<Mesh> coarseMesh,
                  const std::shared_ptr<Mesh> fineMesh);

    //! Builds the hierarchy from coarse to fine mesh
    void build();

    //! Returns the mesh hierarchy
    std::shared_ptr<MeshHierarchyTable> get() const {
        return m_hierarchy;
    }

private:
    const std::shared_ptr<Mesh> m_coarseMesh;
    const std::shared_ptr<Mesh> m_fineMesh;

    std::shared_ptr<MeshHierarchyTable> m_hierarchy;
    std::unique_ptr<PointLocator> m_pLocCoarse;
};

//! Builds mesh hierarchy across multiple-levels
//! using the hierarchies given between successive levels
std::shared_ptr<MeshHierarchyTable> buildMLMeshHierarchy
(std::vector<std::shared_ptr<MeshHierarchyTable>> hierarchies);

//! For a given element in the coarsest mesh,
//! builds mesh hierarchy across multiple-levels
//! using the hierarchies given between successive levels
std::set<int> buildElemMLMeshHierarchy
(int coarseElId,
 std::vector<std::shared_ptr<MeshHierarchyTable>> hierarchies);

}

#endif // MYMFEM_MESH_HIERARCHY_HPP
