#include "nested_hierarchy.hpp"

using namespace mfem;


void mymfem::NestedMeshHierarchy
:: buildHierarchicalTranformations() const
{
    int numMeshes = m_meshes.size();

    m_hierarchicalTransformations.resize(numMeshes-1);

    for (int i=0; i<numMeshes-1; i++)
    {
        assert(m_meshes[i]->GetNE() <= m_meshes[i+1]->GetNE());

        // initialise hierarchical transformation
        m_hierarchicalTransformations[i]
                = std::make_shared<HierarchicalMeshTransformationTable>
                (m_meshes[i]->GetNE());

        buildTranformationBetweenSuccessiveLevels(i);
    }
}

void mymfem::NestedMeshHierarchy
:: buildTranformationBetweenSuccessiveLevels
(int id) const
{
    auto coarseMesh = m_meshes[id];
    auto fineMesh = m_meshes[id+1];

    std::unique_ptr<PointLocator> pointLocatorCoarse
            = std::make_unique<PointLocator>(m_meshes[id].get());

    // number of elements in the fine mesh
    int numElsFine = fineMesh->GetNE();

    int elIdCoarse;
    IntegrationPoint dip; // dummy

    int initElIdCoarse = 0;
    Vector xC(fineMesh->Dimension());

    for (int i=0; i<numElsFine; i++)
    {
        getElementCenter(fineMesh, i, xC);
        std::tie(elIdCoarse, dip)
                = (*pointLocatorCoarse)(xC, initElIdCoarse);

        m_hierarchicalTransformations[id]->add(elIdCoarse, i);
    }
}


// Builds mesh hierarchy across multiple-levels
// using the hierarchies given between successive levels
std::shared_ptr<HierarchicalMeshTransformationTable>
mymfem::buildMultiLevelMeshHierarchy
(mymfem::HierarchicalMeshTransformations& hierMeshTrans)
{
    int numCoarseEls = hierMeshTrans[0]->getNumParents();
    auto multiLevelHierMeshTrans = std::make_shared
            <HierarchicalMeshTransformationTable>(numCoarseEls);

    for (int i=0; i<numCoarseEls; i++) {
        (*multiLevelHierMeshTrans)(i)
                = buildElemMultiLevelMeshHierarchy(i, hierMeshTrans);
    }

    return multiLevelHierMeshTrans;
}

// For a given element in the coarsest mesh,
// builds mesh hierarchy across multiple-levels
// using the hierarchies given between successive levels
std::set<int> mymfem::buildElemMultiLevelMeshHierarchy
(int coarseElId, mymfem::HierarchicalMeshTransformations& hierMeshTrans)
{
    // if only one mesh hierarchy is given,
    // return its children directly
    if (hierMeshTrans.size() == 1) {
        return (*hierMeshTrans[0])(coarseElId);
    }

    // else gather all children at the finest level
    // by recursively calling this function
    std::set<int> allChildren;

    // crete a sub-hierarchy without the first entry from hierarchy
    HierarchicalMeshTransformations
            subHierMeshTrans(hierMeshTrans.begin() + 1,
                             hierMeshTrans.end());

    auto coarseBuf = (*hierMeshTrans[0])(coarseElId);
    for (auto it=coarseBuf.begin(); it!= coarseBuf.end(); it++) {
//        std::cout << "Sub-element " << *it
//                  << " of element " << coarseElId << std::endl;
        auto children
                = buildElemMultiLevelMeshHierarchy(*it, subHierMeshTrans);
        allChildren.insert(children.begin(), children.end());
    }

    return allChildren;
}

// Builds mesh hierarchy across multiple-levels
// using the hierarchical transformation of current mesh
// and the next hierarchy
std::shared_ptr<HierarchicalMeshTransformationTable>
mymfem::buildMultiLevelMeshHierarchy
(std::shared_ptr<HierarchicalMeshTransformationTable>& curHierarchy,
 std::shared_ptr<HierarchicalMeshTransformationTable>& nextHierarchy)
{
    int numCurEls = curHierarchy->getNumParents();
    auto multiLevelHierarchy = std::make_shared
            <HierarchicalMeshTransformationTable>(numCurEls);

    for (int i=0; i<numCurEls; i++)
    {
        auto curBuf = (*curHierarchy)(i);
        for (auto it=curBuf.begin(); it!= curBuf.end(); it++) {
            auto children = (*nextHierarchy)(*it);
            ((*multiLevelHierarchy)(i)).insert(children.begin(),
                                               children.end());
        }
    }

    return multiLevelHierarchy;
}

// End of file
