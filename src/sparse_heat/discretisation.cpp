#include "discretisation.hpp"

#include "assembly.hpp"
#include "utilities.hpp"
#include "../heat/coefficients.hpp"

using namespace mfem;


sparseHeat::LsqSparseXtFem
:: LsqSparseXtFem (const nlohmann::json& config,
                   std::shared_ptr<heat::TestCases>& testCase,
                   const int numLevels,
                   const int minTemporalLevel)
    : m_config (config),
      m_testCase(testCase),
      m_numLevels (numLevels),
      m_minTemporalLevel (minTemporalLevel)
{
    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config, "deg", m_deg, 1);
    if (m_deg > 1) {
        std::cerr << "The sparse implementation has only been implemented "
                     "for polynomial degree 1." << std::endl;
        abort();
    }

    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(m_config, "end_time",
                                        m_endTime, 1);

    m_maxTemporalLevel = m_minTemporalLevel + m_numLevels - 1;
}

sparseHeat::LsqSparseXtFem
:: LsqSparseXtFem
(const nlohmann::json & config,
 std::shared_ptr<heat::TestCases> & testCase,
 const int numLevels,
 const int minTemporalLevel,
 std::shared_ptr<mymfem::NestedMeshHierarchy> & spatialMeshHierachy)
    : LsqSparseXtFem(config, testCase, numLevels, minTemporalLevel)
{
    setNestedFEHierarchyAndSpatialBoundaryDofs(spatialMeshHierachy);
}

sparseHeat::LsqSparseXtFem
:: ~LsqSparseXtFem()
{
    if (m_systemBlock11) {
        delete m_systemBlock11;
    }

    if (m_systemBlock12) {
        delete m_systemBlock12;
    }

    if (m_systemBlock21) {
        delete m_systemBlock21;
    }

    if (m_systemBlock22) {
        delete m_systemBlock22;
    }

    if (m_systemMatrix) {
        delete m_systemMatrix;
    }
}

void sparseHeat::LsqSparseXtFem
:: setNestedFEHierarchyAndSpatialBoundaryDofs
(std::shared_ptr<mymfem::NestedMeshHierarchy>& spatialNestedMeshHierarchy)
{
    setNestedFEHierarchy(spatialNestedMeshHierarchy);
    setSpatialBoundaryDofs();
}

void sparseHeat::LsqSparseXtFem
:: setNestedFEHierarchy
(std::shared_ptr<mymfem::NestedMeshHierarchy>& spatialNestedMeshHierarchy)
{
    setNestedFEHierarchyForTemperature(spatialNestedMeshHierarchy);
    setNestedFEHierarchyForHeatFlux(spatialNestedMeshHierarchy);
}

void sparseHeat::LsqSparseXtFem
:: setNestedFEHierarchyForTemperature
(std::shared_ptr<mymfem::NestedMeshHierarchy>& spatialNestedMeshHierarchy)
{
    auto meshes = spatialNestedMeshHierarchy->getMeshes();
    assert((int) meshes.size() == m_numLevels);

    m_spatialNestedFEHierarchyTemperature
            = std::make_shared<mymfem::NestedFEHierarchy>
            (spatialNestedMeshHierarchy);

    m_xDim = meshes[0]->Dimension();
    FiniteElementCollection *xH1Coll
            = new H1_FECollection(m_deg, m_xDim,
                                  BasisType::GaussLobatto);
    m_xFecs.Append(xH1Coll);

    for (int i=0; i<m_numLevels; i++)
    {
        auto xH1Space = std::make_shared<FiniteElementSpace>
                (meshes[i].get(), xH1Coll);
        m_spatialNestedFEHierarchyTemperature->addFESpace(xH1Space);
    }
}

void sparseHeat::LsqSparseXtFem
:: setSpatialBoundaryDofs()
{
    auto xMeshes
            = m_spatialNestedFEHierarchyTemperature->getMeshes();
    auto xFeSpaces
            = m_spatialNestedFEHierarchyTemperature->getFESpaces();

    auto tHierarchicalFESpacesSizes
            = evalTemporalBlockSizes(m_minTemporalLevel,
                                     m_maxTemporalLevel);
    auto xFeSpacesSizes
            = m_spatialNestedFEHierarchyTemperature->getNumDims();
    auto sizeBuf = evalSpaceTimeBlockSizes(xFeSpacesSizes,
                                           tHierarchicalFESpacesSizes);
    auto offsetsBuf = evalBlockOffsets(sizeBuf);

    for (int m=0; m<m_numLevels; m++)
    {
        int i = getSpatialIndex(m);
        int shift = offsetsBuf[m];
        if (xMeshes[i]->bdr_attributes.Size())
        {
            Array<int> xEssBdrMarker(xMeshes[i]->bdr_attributes.Max());
            m_testCase->setBdryDirichlet(xEssBdrMarker);

            Array<int> xEssDofs;
            xFeSpaces[i]->GetEssentialTrueDofs(xEssBdrMarker, xEssDofs);

            int tNumDofs = tHierarchicalFESpacesSizes[m];
            int xNumDofs = xFeSpacesSizes[i];

            int xNumEssDofs = xEssDofs.Size();
            Array<int> myEssDofs(xNumEssDofs*tNumDofs);
            for (int j=0; j<tNumDofs; j++) {
                for (int k=0; k<xNumEssDofs; k++) {
                    myEssDofs[k + j*xNumEssDofs]
                            = xEssDofs[k] + j*xNumDofs + shift;
                }
            }
            m_essentialDofs.Append(myEssDofs);
//            myEssDofs.Print();
//            std::cout << "\n\n";
        }
    }
}

void sparseHeat::LsqSparseXtFem
:: assembleSystemSubMatrices()
{
    assembleTemporalInitial();
    assembleTemporalMass();
    assembleTemporalStiffness();
    assembleTemporalGradient();

    assembleSpatialMassForTemperature();
    assembleSpatialMassForHeatFlux();
    assembleSpatialStiffnessForTemperature();
    assembleSpatialStiffnessForHeatFlux();
    assembleSpatialGradient();
    assembleSpatialDivergence();
}

void sparseHeat::LsqSparseXtFem
:: assembleTemporalInitial()
{
    auto temporalInitialAssembler
            = std::make_unique<sparseHeat::TemporalInitialMatrixAssembler>
            (m_endTime, m_minTemporalLevel, m_maxTemporalLevel);
    temporalInitialAssembler->assemble();
    m_temporalInitial = temporalInitialAssembler->getBlockMatrix();
}

void sparseHeat::LsqSparseXtFem
:: assembleTemporalMass()
{
    auto temporalMassAssembler
            = std::make_unique<sparseHeat::TemporalMassMatrixAssembler>
            (m_endTime, m_minTemporalLevel, m_maxTemporalLevel);
    temporalMassAssembler->assemble();
    m_temporalMass = temporalMassAssembler->getBlockMatrix();
}

void sparseHeat::LsqSparseXtFem
:: assembleTemporalStiffness()
{
    auto temporalStiffnessAssembler
            = std::make_unique<sparseHeat::TemporalStiffnessMatrixAssembler>
            (m_endTime, m_minTemporalLevel, m_maxTemporalLevel);
    temporalStiffnessAssembler->assemble();
    m_temporalStiffness = temporalStiffnessAssembler->getBlockMatrix();
}

void sparseHeat::LsqSparseXtFem
:: assembleTemporalGradient()
{
    auto temporalGradientAssembler
            = std::make_unique<sparseHeat::TemporalGradientMatrixAssembler>
            (m_endTime, m_minTemporalLevel, m_maxTemporalLevel);
    temporalGradientAssembler->assemble();
    m_temporalGradient = temporalGradientAssembler->getBlockMatrix();
}

void sparseHeat::LsqSparseXtFem
:: assembleSpatialMassForTemperature()
{
    auto spatialMassBilinearForm
            = std::make_unique<mymfem::BlockBilinearForm>
            (m_spatialNestedFEHierarchyTemperature);
    std::shared_ptr<mymfem::BlockBilinearFormIntegrator>
            spatialMassIntegator
            = std::make_shared<sparseHeat::SpatialMassIntegrator>();
    spatialMassBilinearForm->addDomainIntegrator(spatialMassIntegator);
    spatialMassBilinearForm->assemble();
    m_spatialMass1 = spatialMassBilinearForm->getBlockMatrix();
}

void sparseHeat::LsqSparseXtFem
:: assembleSpatialStiffnessForTemperature()
{
    heat::MediumTensorCoeff mediumCoeff(m_testCase);

    auto spatialStiffnessBilinearForm
            = std::make_unique<mymfem::BlockBilinearForm>
            (m_spatialNestedFEHierarchyTemperature);
    std::shared_ptr<mymfem::BlockBilinearFormIntegrator>
            spatialStiffnessIntegator
            = std::make_shared<sparseHeat::SpatialStiffnessIntegrator>
            (mediumCoeff);
    spatialStiffnessBilinearForm->addDomainIntegrator
            (spatialStiffnessIntegator);
    spatialStiffnessBilinearForm->assemble();
    m_spatialStiffness1 = spatialStiffnessBilinearForm->getBlockMatrix();
}

void sparseHeat::LsqSparseXtFem
:: buildSystemMatrix()
{
    buildSystemBlocks();

    Array<int> blockOffsets(3);
    blockOffsets[0] = 0;
    blockOffsets[1] = m_systemBlock11->NumRows();
    blockOffsets[2] = m_systemBlock21->NumRows();
    blockOffsets.PartialSum();

    BlockMatrix buf(blockOffsets);
    buf.SetBlock(0, 0, m_systemBlock11);
    buf.SetBlock(0, 1, m_systemBlock12);
    buf.SetBlock(1, 0, m_systemBlock21);
    buf.SetBlock(1, 1, m_systemBlock22);

    m_systemMatrix = buf.CreateMonolithic();
    applyBCs(*m_systemMatrix);
}

void sparseHeat::LsqSparseXtFem
:: buildSystemBlocks()
{
    buildSystemBlock11();
    buildSystemBlock22();
    buildSystemBlock12();
    buildSystemBlock21();
}

void sparseHeat::LsqSparseXtFem
:: buildSystemBlock11()
{
    auto temporalBlockSizes
            = evalTemporalBlockSizes(m_minTemporalLevel,
                                     m_maxTemporalLevel);
    auto spatialBlockSizes
            = m_spatialNestedFEHierarchyTemperature->getNumDims();

    auto blockSizes = evalSpaceTimeBlockSizes(spatialBlockSizes,
                                              temporalBlockSizes);
    auto blockOffsets = evalBlockOffsets(blockSizes);

    // tmp1 = (K^t_{\ell1, \ell2} + E^t_{\ell1, \ell2})
    //            \otimes M^{x; (1,1)}_{L-\ell1, L-\ell2}
    auto buf = std::make_unique<BlockMatrix>(blockOffsets);
    for (int i=0; i<m_numLevels; i++)
    {
        int ii = getSpatialIndex(i);
        for (int j=0; j<m_numLevels; j++)
        {
            int jj = getSpatialIndex(j);
            auto tmp = Add(m_temporalInitial->GetBlock(i, j),
                           m_temporalStiffness->GetBlock(i, j));
            auto mat = OuterProduct(*tmp, m_spatialMass1->GetBlock(ii, jj));
            buf->SetBlock(i, j, mat);
            delete tmp;
        }
    }
    buf->owns_blocks = true;
    auto tmp1 = buf->CreateMonolithic();

    // tmp2 = M^t_{\ell1, \ell2} \otimes K^{x;(1,1)}_{L-\ell1, L-\ell2}
    buf.reset(new BlockMatrix(blockOffsets));
    for (int i=0; i<m_numLevels; i++)
    {
        int ii = getSpatialIndex(i);
        for (int j=0; j<m_numLevels; j++)
        {
            int jj = getSpatialIndex(j);
            auto mat = OuterProduct(m_temporalMass->GetBlock(i, j),
                                    m_spatialStiffness1->GetBlock(ii, jj));
            buf->SetBlock(i, j, mat);
        }
    }
    buf->owns_blocks = true;
    auto tmp2 = buf->CreateMonolithic();

    m_systemBlock11 = Add(*tmp1, *tmp2);

    delete tmp1;
    delete tmp2;
}

void sparseHeat::LsqSparseXtFem
:: buildSystemBlock22()
{
    auto temporalBlockSizes
            = evalTemporalBlockSizes(m_minTemporalLevel,
                                     m_maxTemporalLevel);
    auto spatialBlockSizes
            = m_spatialNestedFEHierarchyHeatFlux->getNumDims();

    auto blockSizes = evalSpaceTimeBlockSizes(spatialBlockSizes,
                                              temporalBlockSizes);
    auto blockOffsets = evalBlockOffsets(blockSizes);

    auto buf = std::make_unique<BlockMatrix>(blockOffsets);
    // block22 = M^t_{\ell1, \ell2}
    // \otimes (M^{x; (2,2)}_{L-\ell1, L-\ell2}
    //        + K^{x; (2,2)}_{L-\ell1, L-\ell2})
    for (int i=0; i<m_numLevels; i++)
    {
        int ii = getSpatialIndex(i);
        for (int j=0; j<m_numLevels; j++)
        {
            int jj = getSpatialIndex(j);
            auto tmp = Add(m_spatialMass2->GetBlock(ii, jj),
                           m_spatialStiffness2->GetBlock(ii, jj));
            auto mat = OuterProduct(m_temporalMass->GetBlock(i, j), *tmp);
            buf->SetBlock(i, j, mat);
            delete tmp;
        }
    }
    buf->owns_blocks = true;
    m_systemBlock22 = buf->CreateMonolithic();
}

void sparseHeat::LsqSparseXtFem
:: buildSystemBlock12()
{
    auto temporalBlockSizes
            = evalTemporalBlockSizes(m_minTemporalLevel,
                                     m_maxTemporalLevel);
    auto spatialRowBlockSizes
            = m_spatialNestedFEHierarchyTemperature->getNumDims();
    auto spatialColBlockSizes
            = m_spatialNestedFEHierarchyHeatFlux->getNumDims();

    auto rowBlockSizes = evalSpaceTimeBlockSizes(spatialRowBlockSizes,
                                                 temporalBlockSizes);
    auto colBlockSizes = evalSpaceTimeBlockSizes(spatialColBlockSizes,
                                                 temporalBlockSizes);

    auto rowBlockOffsets = evalBlockOffsets(rowBlockSizes);
    auto colBlockOffsets = evalBlockOffsets(colBlockSizes);

    auto buf = std::make_unique<BlockMatrix>(rowBlockOffsets,
                                             colBlockOffsets);
    // tmp1 = (C^t_{\ell1, \ell2})^{\top}
    //           \otimes B^{x; (1,2)}_{L-\ell1, L-\ell2})
    for (int i=0; i<m_numLevels; i++)
    {
        int ii = getSpatialIndex(i);
        for (int j=0; j<m_numLevels; j++)
        {
            int jj = getSpatialIndex(j);
            auto tmp = Transpose(m_temporalGradient->GetBlock(j, i));
            auto mat = OuterProduct
                    (*tmp, m_spatialDivergence->GetBlock(ii, jj));
//            std::cout << i << "\t"
//                      << j << "\t"
//                      << mat->NumRows() << "\t"
//                      << mat->NumCols() << "\t"
//                      << rowBlockSizes[i] << "\t"
//                      << colBlockSizes[j] << std::endl;
            buf->SetBlock(i, j, mat);
            delete tmp;
        }
    }
    buf->owns_blocks = true;
    auto tmp1 = buf->CreateMonolithic();

    // tmp2 = (M^t_{\ell1, \ell2})^{\top}
    //           \otimes C^{x; (1,2)}_{L-\ell1, L-\ell2}
    buf.reset(new BlockMatrix(rowBlockOffsets, colBlockOffsets));
    for (int i=0; i<m_numLevels; i++)
    {
        int ii = getSpatialIndex(i);
        for (int j=0; j<m_numLevels; j++)
        {
            int jj = getSpatialIndex(j);
            auto tmp = Transpose(m_spatialGradient->GetBlock(jj, ii));
            auto mat = OuterProduct(m_temporalMass->GetBlock(i, j), *tmp);
            buf->SetBlock(i, j, mat);
            delete tmp;
        }
    }
    buf->owns_blocks = true;
    auto tmp2 = buf->CreateMonolithic();

    m_systemBlock12 = Add(-1., *tmp1, -1., *tmp2);

    delete tmp1;
    delete tmp2;
}

void sparseHeat::LsqSparseXtFem
:: buildSystemBlock21()
{
    m_systemBlock21 = Transpose(*m_systemBlock12);
}

int sparseHeat::LsqSparseXtFem
:: getSpatialIndex(int i) const {
    return m_numLevels - i - 1;
}

void sparseHeat::LsqSparseXtFem
:: assembleRhs(std::shared_ptr<BlockVector>& B) const
{
    Vector &B0 = B->GetBlock(0);
    assembleICs(B0);

    assembleSource(B);

    applyBCs(*B);
}

void sparseHeat::LsqSparseXtFem
:: assembleICs(Vector& b) const
{
    auto xFeSpaces
            = m_spatialNestedFEHierarchyTemperature->getFESpaces();
    int coarsestTemporalIndex = 0;
    int finestSpatialIndex = getSpatialIndex(coarsestTemporalIndex);

    LinearForm icsForm(xFeSpaces[finestSpatialIndex].get());
    heat::InitialTemperatureCoeff u0(m_testCase);
    icsForm.AddDomainIntegrator
            (new DomainLFIntegrator(u0));
    icsForm.Assemble();
    Vector buf(icsForm.GetData(),
               xFeSpaces[finestSpatialIndex]->GetTrueVSize());

    b.SetVector(buf, 0);
}

void sparseHeat::LsqSparseXtFem
:: assembleSource(std::shared_ptr<BlockVector>& B) const
{
    assembleSourceWithTemporalGradientOfTemperatureBasis(B->GetBlock(0));
    assembleSourceWithSpatialDivergenceOfHeatFluxBasis(B->GetBlock(1));
}

void sparseHeat::LsqSparseXtFem
:: assembleSourceWithTemporalGradientOfTemperatureBasis(Vector& b) const
{
    auto spatialFes = m_spatialNestedFEHierarchyTemperature->getFESpaces();

    // auxiliary variables
    Vector buf;
    int shift = 0;

    // temporal integration rule
    int order = 4;
    const IntegrationRule *ir
            = &IntRules.Get(Geometry::SEGMENT, order);

    // coarsest temporal + finest spatial level
    // standard FE basis at the coarsest temporal level
    {
        int numTemporalMeshElements
                = static_cast<int>(std::pow(2, m_minTemporalLevel));
        double ht = m_endTime / numTemporalMeshElements;

        int spatialId = getSpatialIndex(0);
        auto curSpatialFes = spatialFes[spatialId];
        int curSpatialFesDim = curSpatialFes->GetTrueVSize();

        buf.SetSize(curSpatialFesDim);
        for (int n=0; n<numTemporalMeshElements; n++)
        {
            double tLeft = n*ht;
            double tRight = (n+1)*ht;

            Vector tmp1(b.GetData() + n*curSpatialFesDim + shift,
                        curSpatialFesDim);
            Vector tmp2(b.GetData() + (n+1)*curSpatialFesDim + shift,
                        curSpatialFesDim);

            for (int i=0; i<ir->GetNPoints(); i++)
            {
                const IntegrationPoint &ip = ir->IntPoint(i);
                double t = affineTransform(tLeft, tRight, ip.x);

                assembleSourceWithTemporalGradientOfTemperatureBasisAtGivenTime
                        (t, curSpatialFes, buf);

                // tmp1 += weight * temporalGradient1 * ht * ft
                // temporalGradient1 * ht = -1
                tmp1.Add(-ip.weight, buf);

                // tmp2 += weight * temporalGradient2 * ht * ft
                // temporalGradient2 * ht = +1
                tmp2.Add(+ip.weight, buf);

//                std::cout << n << "\t"
//                          << tLeft << "\t"
//                          << tRight << "\t"
//                          << i << "\t"
//                          << t << std::endl;
            }
        }
        shift += (numTemporalMeshElements+1)*curSpatialFesDim;
    }

    // hierarchical temporal basis for the remaining levels
    for (int m=1; m<m_numLevels; m++)
    {
        int numTemporalMeshElements
                = static_cast<int>(std::pow(2, m_minTemporalLevel + m));
        double ht = m_endTime / numTemporalMeshElements;

        int spatialId = getSpatialIndex(m);
        auto curSpatialFes = spatialFes[spatialId];
        int curSpatialFesDim = curSpatialFes->GetTrueVSize();

        buf.SetSize(curSpatialFesDim);
        for (int n=0; n<numTemporalMeshElements; n++)
        {
            double tLeft = n*ht;
            double tRight = (n+1)*ht;

            int nn = n/2;
            Vector tmp(b.GetData() + nn*curSpatialFesDim + shift, curSpatialFesDim);

            for (int i=0; i<ir->GetNPoints(); i++)
            {
                const IntegrationPoint &ip = ir->IntPoint(i);
                double t = affineTransform(tLeft, tRight, ip.x);

                assembleSourceWithTemporalGradientOfTemperatureBasisAtGivenTime
                        (t, curSpatialFes, buf);

                if (n%2 == 0) {
                    // tmp += weight * temporalGradient * ht * f(t)
                    // temporalGradient * h = +1
                    tmp.Add(+ip.weight, buf);
                }
                else {
                    // tmp += weight * temporalGradient * ht * ft
                    // temporalGradient * h = -1
                    tmp.Add(-ip.weight, buf);
                }
            }
        }
        shift += (numTemporalMeshElements/2)*curSpatialFesDim;
    }
}

void sparseHeat::LsqSparseXtFem
:: assembleSourceWithTemporalGradientOfTemperatureBasisAtGivenTime
(double t,
 const std::shared_ptr<FiniteElementSpace>& spatialFes,
 Vector& b) const
{
    LinearForm sourceForm(spatialFes.get());
    heat::SourceCoeff sourceCoeff(m_testCase);
    sourceCoeff.SetTime(t);
    sourceForm.AddDomainIntegrator
            (new DomainLFIntegrator(sourceCoeff, 2, 4));
    sourceForm.Assemble();
    b = sourceForm.GetData();
}

void sparseHeat::LsqSparseXtFem
:: assembleSourceWithSpatialDivergenceOfHeatFluxBasis(Vector& b) const
{
    auto spatialFes = m_spatialNestedFEHierarchyHeatFlux->getFESpaces();

    // auxiliary variables
    Vector buf;
    int shift = 0;

    // temporal integration rule
    int order = 4;
    const IntegrationRule *ir
            = &IntRules.Get(Geometry::SEGMENT, order);

    // coarsest temporal + finest spatial level
    // standard FE basis at the coarsest temporal level
    {
        int numTemporalMeshElements
                = static_cast<int>(std::pow(2, m_minTemporalLevel));
        double ht = m_endTime / numTemporalMeshElements;

        int spatialId = getSpatialIndex(0);
        auto curSpatialFes = spatialFes[spatialId];
        int curSpatialFesDim = curSpatialFes->GetTrueVSize();

        buf.SetSize(curSpatialFesDim);
        for (int n=0; n<numTemporalMeshElements; n++)
        {
            double tLeft = n*ht;
            double tRight = (n+1)*ht;

            Vector tmp1(b.GetData() + n*curSpatialFesDim + shift,
                        curSpatialFesDim);
            Vector tmp2(b.GetData() + (n+1)*curSpatialFesDim + shift,
                        curSpatialFesDim);

            for (int i=0; i<ir->GetNPoints(); i++)
            {
                const IntegrationPoint &ip = ir->IntPoint(i);
                double t = affineTransform(tLeft, tRight, ip.x);

                assembleSourceWithSpatialDivergenceOfHeatFluxBasisAtGivenTime
                        (t, curSpatialFes, buf);

                // tmp1 += weight * temporalBasisVal1 * ht * ft
                tmp1.Add(ip.weight * evalRightHalfOfHatBasis(t, tLeft, ht) * ht, buf);

                // tmp2 += weight * temporalBasisVal2 * ht * ft
                tmp2.Add(+ip.weight * evalLeftHalfOfHatBasis(t, tLeft, ht) * ht, buf);
            }
        }
        shift += (numTemporalMeshElements+1)*curSpatialFesDim;
    }

    // hierarchical temporal basis for the remaining levels
    for (int m=1; m<m_numLevels; m++)
    {
        int numTemporalMeshElements
                = static_cast<int>(std::pow(2, m_minTemporalLevel + m));
        double ht = m_endTime / numTemporalMeshElements;

        int spatialId = getSpatialIndex(m);
        auto curSpatialFes = spatialFes[spatialId];
        int curSpatialFesDim = curSpatialFes->GetTrueVSize();

        buf.SetSize(curSpatialFesDim);
        for (int n=0; n<numTemporalMeshElements; n++)
        {
            double tLeft = n*ht;
            double tRight = (n+1)*ht;

            int nn = n/2;
            Vector tmp(b.GetData() + nn*curSpatialFesDim + shift, curSpatialFesDim);

            for (int i=0; i<ir->GetNPoints(); i++)
            {
                const IntegrationPoint &ip = ir->IntPoint(i);
                double t = affineTransform(tLeft, tRight, ip.x);

                assembleSourceWithSpatialDivergenceOfHeatFluxBasisAtGivenTime
                        (t, curSpatialFes, buf);

                if (n%2 == 0) {
                    // tmp += weight * temporalBasisVal * ht * f(t)
                    tmp.Add(ip.weight * evalLeftHalfOfHatBasis(t, tLeft, ht) * ht, buf);
                }
                else {
                    // tmp += weight * temporalGradient * ht * ft
                    tmp.Add(ip.weight * evalRightHalfOfHatBasis(t, tLeft, ht) * ht, buf);
                }
            }
        }
        shift += (numTemporalMeshElements/2)*curSpatialFesDim;
    }
}

void sparseHeat::LsqSparseXtFem
:: applyBCs(SparseMatrix &A) const
{
    for (int k=0; k<m_essentialDofs.Size(); k++) {
        A.EliminateRowCol(m_essentialDofs[k]);
    }
}

void sparseHeat::LsqSparseXtFem
:: applyBCs(BlockVector &B) const
{
    auto& b = B.GetBlock(0);
    b.SetSubVector(m_essentialDofs, 0.);
}

// End of file
