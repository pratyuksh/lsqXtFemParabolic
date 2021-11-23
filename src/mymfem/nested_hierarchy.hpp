#ifndef MYMFEM_NESTED_HIERARCHY_HPP
#define MYMFEM_NESTED_HIERARCHY_HPP

#include "mfem.hpp"
using namespace mfem;

#include <vector>
#include <set>
#include <memory>

#include "utilities.hpp"


typedef struct HierarchicalMeshTransformationTable
{
public:
    //! Constructor with number of coarse mesh elements
    HierarchicalMeshTransformationTable(int numElsCoarse)
        : m_numRows(numElsCoarse) {
        m_data.resize(m_numRows);
    }

    //! Adds an element in the fine mesh to the hierarchy of an element
    //! in the coarse mesh
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
} HierarchicalMeshTransformationTable;


namespace mymfem
{
typedef std::vector<std::shared_ptr<Mesh>> NestedMeshes;
typedef std::vector<std::shared_ptr<FiniteElementSpace>> NestedFESpaces;
typedef std::vector<std::shared_ptr<HierarchicalMeshTransformationTable>>
HierarchicalMeshTransformations;

/**
 * @brief Wrapper for a sequence of nested meshes.
 * Also provides routines to build the hierarchical transformations
 * between successive meshes
 */
struct NestedMeshHierarchy {
public:
    //! Default constructor
    NestedMeshHierarchy() {}

    //! Adds mesh to the hierarchy
    void addMesh(std::shared_ptr<Mesh>& mesh) {
        m_meshes.push_back(mesh);
    }

    void finalize() {
        buildHierarchicalTranformations();
    }

    //! Builds the hierarchy transformations
    //! between all successive meshes
    void buildHierarchicalTranformations() const;

    //! Returns the meshes
    NestedMeshes& getMeshes() {
      return m_meshes;
    }

    NestedMeshes getMeshes() const {
      return m_meshes;
    }

    //! Returns the number of meshes in the hierarchy
    int getNumMeshes() const {
        return std::move(m_meshes.size());
    }

    //! Returns the array of hierarchical mesh transformations
    HierarchicalMeshTransformations getTransformations() const {
        return m_hierarchicalTransformations;
    }

private:
    //! Builds the hierarchy tranformation between two successive meshes
    void buildTranformationBetweenSuccessiveLevels (int id) const;

private:
    NestedMeshes m_meshes;
    mutable HierarchicalMeshTransformations m_hierarchicalTransformations;
};

//! Builds mesh hierarchy across multiple-levels
//! using the hierarchies given between successive levels
std::shared_ptr<HierarchicalMeshTransformationTable>
buildMultiLevelMeshHierarchy
(HierarchicalMeshTransformations& hierMeshTrans);

//! For a given element in the coarsest mesh,
//! builds mesh hierarchy across multiple-levels
//! using the hierarchies given between successive levels
std::set<int> buildElemMultiLevelMeshHierarchy
(int coarseElId, HierarchicalMeshTransformations& hierMeshTrans);

//! Builds mesh hierarchy across multiple-levels
//! using the hierarchical transformation of current mesh
//! and the next hierarchy
std::shared_ptr<HierarchicalMeshTransformationTable>
buildMultiLevelMeshHierarchy
(std::shared_ptr<HierarchicalMeshTransformationTable>& curHierarchy,
 std::shared_ptr<HierarchicalMeshTransformationTable>& nextHierarchy);


/**
 * @brief Wrapper for a sequence of nested finite element spaces.
 * Contains both the nested mesh hierarchy
 * and the nested finite element spaces.
 */
struct NestedFEHierarchy {
public:
    //! Default constructor
    NestedFEHierarchy(std::shared_ptr<NestedMeshHierarchy>&
                      nestedMeshHierarchy)
      : m_nestedMeshHierarchy(nestedMeshHierarchy) {}

    //! Adds FE space to the hierarchy
    void addFESpace(std::shared_ptr<FiniteElementSpace>& fes) {
        m_feSpaces.push_back(fes);
        
        int n = m_feSpaces.size();
        auto meshes = m_nestedMeshHierarchy->getMeshes();
        assert(m_feSpaces[n-1]->GetMesh() == meshes[n-1].get());
    }

    //! Returns the number of levels in the hierarchy
    int getNumLevels() const {
        return std::move(m_nestedMeshHierarchy->getNumMeshes());
    }

    //! Returns the nested FE spaces
    NestedFESpaces getFESpaces() const {
        return m_feSpaces;
    }

    //! Returns the nested meshes
    NestedMeshes getMeshes() {
        return m_nestedMeshHierarchy->getMeshes();
    }

    //! Returns the hierarchical mesh transformations
    HierarchicalMeshTransformations getHierarchicalMeshTransformations() {
        return m_nestedMeshHierarchy->getTransformations();
    }

    //! Returns the number of dimensions (size) of the FE spaces
    Array<int> getNumDims()
    {
        int numLevels = this->getNumLevels();
        Array<int> numDims(numLevels);
        for (int i=0; i < numLevels; i++) {
            numDims[i] = m_feSpaces[i]->GetTrueVSize();
        }

        return numDims;
    }

private:
    NestedFESpaces m_feSpaces;
    std::shared_ptr<NestedMeshHierarchy> m_nestedMeshHierarchy;
};

}


#endif // MYMFEM_NESTED_HIERARCHY_HPP
