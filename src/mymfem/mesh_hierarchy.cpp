#include "mesh_hierarchy.hpp"

mymfem::MeshHierarchy
:: MeshHierarchy (const std::shared_ptr<Mesh> coarseMesh,
                  const std::shared_ptr<Mesh> fineMesh)
    : m_coarseMesh(coarseMesh), m_fineMesh(fineMesh)
{
    assert(m_coarseMesh->GetNE() <= m_fineMesh->GetNE());

    // initialise mesh hierarchy
    m_hierarchy = std::make_shared<MeshHierarchyTable>
            (m_coarseMesh->GetNE());

    // initialise point locator
    m_pLocCoarse = std::make_unique<PointLocator> (coarseMesh.get());
}

void mymfem::MeshHierarchy
:: build()
{
    // number of elements in the fine mesh
    int numElsFine = m_fineMesh->GetNE();

    int elIdCoarse;
    IntegrationPoint dip; // dummy

    int initElIdCoarse = 0;
    Vector xC(m_fineMesh->Dimension());

    for (int i=0; i<numElsFine; i++)
    {
        getElementCenter(m_fineMesh, i, xC);
        std::tie(elIdCoarse, dip) = (*m_pLocCoarse)(xC, initElIdCoarse);

        m_hierarchy->add(elIdCoarse, i);
    }
}


// Builds mesh hierarchy across multiple-levels
// using the hierarchies given between successive levels
std::shared_ptr<MeshHierarchyTable> mymfem::buildMLMeshHierarchy
(std::vector<std::shared_ptr<MeshHierarchyTable>> hierarchies)
{
    int numCoarseEls = hierarchies[0]->getNumRows();
    auto multiLevelHierarchy
            = std::make_shared<MeshHierarchyTable>(numCoarseEls);

    for (int i=0; i<numCoarseEls; i++){
        (*multiLevelHierarchy)(i)
                = buildElemMLMeshHierarchy(i, hierarchies);
    }

    return multiLevelHierarchy;
}

// For a given element in the coarsest mesh,
// builds mesh hierarchy across multiple-levels
// using the hierarchies given between successive levels
std::set<int> mymfem::buildElemMLMeshHierarchy
(int coarseElId,
 std::vector<std::shared_ptr<MeshHierarchyTable>> hierarchies)
{
    // if only one mesh hierarchy is given,
    // return its children directly
    if (hierarchies.size() == 1) {
        return (*hierarchies[0])(coarseElId);
    }

    // else gather all children at the finest level
    // by recursively calling this function
    std::set<int> allChildren;

    // crete a sub-hierarchy without the first entry from hierarchies
    std::vector<std::shared_ptr<MeshHierarchyTable>>
            subHierarchies(hierarchies.begin() + 1, hierarchies.end());

    auto coarseBuf = (*hierarchies[0])(coarseElId);
    for (auto it=coarseBuf.begin(); it!= coarseBuf.end(); it++) {
        std::cout << "Sub-element " << *it
                  << " of element " << coarseElId << std::endl;
        auto children = buildElemMLMeshHierarchy(*it, subHierarchies);
        allChildren.insert(children.begin(), children.end());
    }

    return allChildren;
}

// End of file
