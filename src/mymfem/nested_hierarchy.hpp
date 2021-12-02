#ifndef MYMFEM_NESTED_HIERARCHY_HPP
#define MYMFEM_NESTED_HIERARCHY_HPP

#include "mfem.hpp"

#include <vector>
#include <set>
#include <memory>

#include "utilities.hpp"


typedef struct HierarchicalMeshTransformationTable
{
public:
    //! Constructor with number of parent mesh elements
    HierarchicalMeshTransformationTable(int numParentEls)
        : m_numParents(numParentEls) {
        m_data.resize(m_numParents);
    }

    //! Adds a child element to the hierarchy of a parent element
    void add (int parentElId, int childElId) {
        m_data[parentElId].insert(childElId);
    }
   
    //! Returns the child indices for the parent element with
    //! index parentElId
    std::set<int>& operator() (int parentElId) {
        return m_data[parentElId];
    }

    //! Returns the number of parent elements
    int getNumParents() const {
        return m_numParents;
    }

    //! Returns the number of children for the parent element with
    //! index parentElId
    int getNumChildren(int parentElId) const {
        return m_data[parentElId].size();
    }

    //! Returns the index of the parent element for a given child element
    int getParentId(int childId)
    {
        for (int i=0; i<m_numParents; i++) {
          if (m_data[i].find(childId) != m_data[i].end()) {
              return i;
          }
        }
#ifndef NDEBUG
        std::cout << "Child element with id " << childId
                  << " not found in the hierarchical "
                     "mesh transformation table!\n";
#endif
        return -1;
    }

    void print() {
        for (int i=0; i<m_numParents; i++) {
            for (auto it=m_data[i].begin(); it!= m_data[i].end(); it++) {
                std::cout << *it << "\t";
            }
            std::cout << "\n";
        }
    }

private:
   int m_numParents; // number of parent mesh elements
   std::vector<std::set<int>> m_data; // stores the mesh hierarchy
} HierarchicalMeshTransformationTable;


namespace mymfem
{
typedef std::vector<std::shared_ptr<mfem::Mesh>>
NestedMeshes;

typedef std::vector<std::shared_ptr<mfem::FiniteElementSpace>>
NestedFESpaces;

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
    void addMesh(std::shared_ptr<mfem::Mesh>& mesh) {
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
    void addFESpace(std::shared_ptr<mfem::FiniteElementSpace>& fes) {
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

    //! Returns the nested mesh hierarchy
    std::shared_ptr<NestedMeshHierarchy> getNestedMeshHierarchy() const {
        return m_nestedMeshHierarchy;
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
    mfem::Array<int> getNumDims();

private:
    NestedFESpaces m_feSpaces;
    std::shared_ptr<NestedMeshHierarchy> m_nestedMeshHierarchy;
};

}


#endif // MYMFEM_NESTED_HIERARCHY_HPP
