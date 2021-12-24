#include "discretisation.hpp"
#include "coefficients.hpp"
#include "assembly.hpp"
#include "../mymfem/utilities.hpp"

#include <fstream>

using namespace mfem;


heat::LsqXtFem
:: LsqXtFem (const nlohmann::json& config,
             std::shared_ptr<heat::TestCases>& testCase)
    : m_config (config),
      m_testCase(testCase)
{
    READ_CONFIG_PARAM_OR_SET_TO_DEFAULT(config, "deg", m_deg, 1);
}

heat::LsqXtFem :: ~LsqXtFem()
{
    reset();
}

// Releases all the memory allocated by this class
void heat::LsqXtFem
:: reset(bool firstPass)
{
    m_firstPassOfAssembleSystemSubMatrices = firstPass;
    m_firstPassOfAssembleSystemBlocks = firstPass;

    resetFeSpaces();
    resetSystemOperators();
    resetSystemBlocks();
    resetSystemSubMatrices();
}

void heat::LsqXtFem
:: resetFeSpaces()
{
    if (m_temporalFeCollection) {
        delete m_temporalFeCollection;
        m_temporalFeCollection = nullptr;
    }
    if (m_temporalFeSpace) {
        delete m_temporalFeSpace;
        m_temporalFeSpace = nullptr;
    }

    for (int i=0; i<m_spatialFeCollections.Size(); i++)
    {
        if (m_spatialFeCollections[i]) {
            delete m_spatialFeCollections[i];
            m_spatialFeCollections[i] = nullptr;
        }
        if (m_spatialFeSpaces[i]) {
            delete m_spatialFeSpaces[i];
            m_spatialFeSpaces[i] = nullptr;
        }
    }
    m_spatialFeCollections.DeleteAll();
}

void heat::LsqXtFem
:: resetSystemOperators()
{
    if (m_systemOperator) {
        delete m_systemOperator;
        m_systemOperator = nullptr;
    }
    if (m_systemMatrix) {
        delete m_systemMatrix;
        m_systemMatrix = nullptr;
    }
}

void heat::LsqXtFem
:: resetSystemBlocks()
{
    resetMediumIndependentSystemBlocks();
    resetMediumDependentSystemBlocks();
}

void heat::LsqXtFem
:: resetMediumIndependentSystemBlocks()
{
    if (m_systemBlock22) {
        clear(m_systemBlock22);
    }

    if (m_materialIndependentSystemBlock11) {
        clear(m_materialIndependentSystemBlock11);
    }
    if (m_materialIndependentSystemBlock12) {
        clear(m_materialIndependentSystemBlock12);
    }
}

void heat::LsqXtFem
:: resetMediumDependentSystemBlocks()
{
    if (m_systemBlock11) {
        clear(m_systemBlock11);
    }
    if (m_systemBlock12) {
        clear(m_systemBlock12);
    }
    if (m_systemBlock21) {
        clear(m_systemBlock21);
    }
}

void heat::LsqXtFem
:: resetSystemSubMatrices()
{
    resetMediumIndependentSystemSubMatrices();
    resetMediumDependentSystemSubMatrices();
}

void heat::LsqXtFem
:: resetMediumIndependentSystemSubMatrices()
{
    if (m_temporalInitial) {
        clear(m_temporalInitial);
    }
    if (m_temporalMass) {
        clear(m_temporalMass);
    }
    if (m_temporalStiffness) {
        clear(m_temporalStiffness);
    }
    if (m_temporalGradient) {
        clear(m_temporalGradient);
    }

    if (m_spatialMass1) {
        clear(m_spatialMass1);
    }
    if (m_spatialMass2) {
        clear(m_spatialMass2);
    }
    if (m_spatialStiffness2) {
        clear(m_spatialStiffness2);
    }
    if (m_spatialDivergence) {
        clear(m_spatialDivergence);
    }
}

void heat::LsqXtFem
:: resetMediumDependentSystemSubMatrices()
{
    if (m_spatialStiffness1) {
        clear(m_spatialStiffness1);
    }
    if (m_spatialGradient) {
        clear(m_spatialGradient);
    }
}

void heat::LsqXtFem
:: resetFeSpacesAndBlockOffsetsAndSpatialBoundaryDofs
(std::shared_ptr<Mesh> &temporalMesh, std::shared_ptr<Mesh> &spatialMesh)
{
    resetFeSpaces();
    setFeSpacesAndBlockOffsetsAndSpatialBoundaryDofs(temporalMesh,
                                                     spatialMesh);
}

void heat::LsqXtFem
:: setFeSpacesAndBlockOffsetsAndSpatialBoundaryDofs
(std::shared_ptr<Mesh> &temporalMesh, std::shared_ptr<Mesh> &spatialMesh)
{
    setFeSpaces(temporalMesh, spatialMesh);
    setBlockOffsets();
    setSpatialBoundaryDofs();
}

void heat::LsqXtFem
:: setFeSpaces(std::shared_ptr<Mesh> &temporalMesh,
               std::shared_ptr<Mesh> &spatialMesh)
{
    // temporal FE space
    m_temporalFeCollection = new H1_FECollection(m_deg, 1,
                                 BasisType::GaussLobatto);
    m_temporalFeSpace = new FiniteElementSpace(temporalMesh.get(),
                                        m_temporalFeCollection);

    // spatial FE spaces
    setSpatialFeSpaceForTemperature(spatialMesh);
    setSpatialFeSpaceForHeatFlux(spatialMesh);
}

void heat::LsqXtFem
:: setSpatialFeSpaceForTemperature(std::shared_ptr<Mesh> &spatialMesh)
{
    m_xDim = spatialMesh->Dimension();
    FiniteElementCollection *xH1Coll
            = new H1_FECollection(m_deg, m_xDim,
                                  BasisType::GaussLobatto);
    auto spatialFeSpaceForTemperature
            = new FiniteElementSpace(spatialMesh.get(), xH1Coll);

    m_spatialFeCollections.Append(xH1Coll);
    m_spatialFeSpaces.Append(spatialFeSpaceForTemperature);
}

void heat::LsqXtFem
:: setBlockOffsets()
{
    int tFeSpaceSize = m_temporalFeSpace->GetTrueVSize();
    int xFeSpaceSizeForTemperature = m_spatialFeSpaces[0]->GetTrueVSize();
    int xFeSpaceSizeForHeatFlux = m_spatialFeSpaces[1]->GetTrueVSize();

    m_blockOffsets.SetSize(3);
    m_blockOffsets[0] = 0;
    m_blockOffsets[1] = xFeSpaceSizeForTemperature*tFeSpaceSize;
    m_blockOffsets[2] = xFeSpaceSizeForHeatFlux*tFeSpaceSize;
    m_blockOffsets.PartialSum();
#ifndef NDEBUG
    std::cout << "\tBlock offsets:" << "\t";
    m_blockOffsets.Print();
#endif
}

void heat::LsqXtFem
:: setSpatialBoundaryDofs()
{
    auto spatialMesh = m_spatialFeSpaces[0]->GetMesh();
    int tFeSpaceSize = m_temporalFeSpace->GetTrueVSize();
    int xFeSpaceSize = m_spatialFeSpaces[0]->GetTrueVSize();

    if (spatialMesh->bdr_attributes.Size())
    {
        m_spatialEssentialBoundaryMarker.SetSize
                (spatialMesh->bdr_attributes.Max());
        m_testCase->setBdryDirichlet(m_spatialEssentialBoundaryMarker);

        m_spatialFeSpaces[0]->GetEssentialTrueDofs
                (m_spatialEssentialBoundaryMarker, m_spatialEssentialDofs);

        int numSpatialEssentialDofs = m_spatialEssentialDofs.Size();
        m_essentialDofs.SetSize(numSpatialEssentialDofs*tFeSpaceSize);
        for (int j=0; j<tFeSpaceSize; j++)
            for (int i=0; i<numSpatialEssentialDofs; i++)
                m_essentialDofs[i + j*numSpatialEssentialDofs]
                        = j*xFeSpaceSize
                        + m_spatialEssentialDofs[i];
    }
}

void heat::LsqXtFem
:: reassembleSystemSubMatrices()
{
    if (m_firstPassOfAssembleSystemSubMatrices) {
        assembleMaterialIndependentSystemSubMatrices();
    }

    resetMediumDependentSystemSubMatrices();
    assembleMaterialDependentSystemSubMatrices();

    if (m_firstPassOfAssembleSystemSubMatrices) {
        m_firstPassOfAssembleSystemSubMatrices = false;
    }
}

void heat::LsqXtFem
:: assembleSystemSubMatrices()
{
    assembleMaterialIndependentSystemSubMatrices();
    assembleMaterialDependentSystemSubMatrices();
}

void heat::LsqXtFem
:: assembleMaterialIndependentSystemSubMatrices()
{
    assembleTemporalInitial();
    assembleTemporalMass();
    assembleTemporalStiffness();
    assembleTemporalGradient();

    assembleSpatialMassForTemperature();
    assembleSpatialMassForHeatFlux();
    assembleSpatialStiffnessForHeatFlux();
    assembleSpatialDivergence();
}

void heat::LsqXtFem
:: assembleMaterialDependentSystemSubMatrices()
{
    assembleSpatialStiffnessForTemperature();
    assembleSpatialGradient();
}

void heat::LsqXtFem
:: assembleTemporalInitial()
{
    Vector tCanonicalBasis(m_temporalFeSpace->GetVSize());
    tCanonicalBasis = 0.0;
    tCanonicalBasis(0) = 1;
    m_temporalInitial = new SparseMatrix(tCanonicalBasis);
}

void heat::LsqXtFem
:: assembleTemporalMass()
{
    BilinearForm *temporalMassForm = new BilinearForm(m_temporalFeSpace);
    temporalMassForm->AddDomainIntegrator(new MassIntegrator);
    temporalMassForm->Assemble();
    temporalMassForm->Finalize();
    m_temporalMass = temporalMassForm->LoseMat();
    delete temporalMassForm;
}

void heat::LsqXtFem
:: assembleTemporalStiffness()
{
    BilinearForm *temporalStiffnessForm = new BilinearForm(m_temporalFeSpace);
    temporalStiffnessForm->AddDomainIntegrator(new DiffusionIntegrator);
    temporalStiffnessForm->Assemble();
    temporalStiffnessForm->Finalize();
    m_temporalStiffness = temporalStiffnessForm->LoseMat();
    delete temporalStiffnessForm;
}

void heat::LsqXtFem
:: assembleTemporalGradient()
{
    BilinearForm *temporalGradientForm
            = new BilinearForm(m_temporalFeSpace);
    temporalGradientForm->AddDomainIntegrator
            (new heat::TemporalGradientIntegrator);
    temporalGradientForm->Assemble();
    temporalGradientForm->Finalize();
    m_temporalGradient = temporalGradientForm->LoseMat();
    delete temporalGradientForm;
}

void heat::LsqXtFem
:: assembleSpatialMassForTemperature()
{
    BilinearForm *spatialMassForm
            = new BilinearForm(m_spatialFeSpaces[0]);
    spatialMassForm->AddDomainIntegrator(new MassIntegrator);
    spatialMassForm->Assemble();
    spatialMassForm->Finalize();
    m_spatialMass1 = spatialMassForm->LoseMat();
    delete spatialMassForm;
}

void heat::LsqXtFem
:: assembleSpatialStiffnessForTemperature()
{
    heat::MediumTensorCoeff mediumCoeff(m_testCase);

    BilinearForm *spatialStiffnessForm
            = new BilinearForm(m_spatialFeSpaces[0]);
    spatialStiffnessForm->AddDomainIntegrator
            (new heat::SpatialStiffnessIntegrator(&mediumCoeff));
    spatialStiffnessForm->Assemble();
    spatialStiffnessForm->Finalize();
    m_spatialStiffness1 = spatialStiffnessForm->LoseMat();
    delete spatialStiffnessForm;
}

// Builds the linear system matrix from the blocks as an MFEM operator
void heat::LsqXtFem
:: rebuildSystemOperator()
{
    resetSystemOperators();

    rebuildSystemBlocks();

    m_systemOperator = new BlockOperator(m_blockOffsets);
    m_systemOperator->SetBlock(0,0, m_systemBlock11);
    m_systemOperator->SetBlock(0,1, m_systemBlock12);
    m_systemOperator->SetBlock(1,0, m_systemBlock21);
    m_systemOperator->SetBlock(1,1, m_systemBlock22);
}

void heat::LsqXtFem
:: buildSystemOperator()
{
    buildSystemBlocks();

    m_systemOperator = new BlockOperator(m_blockOffsets);
    m_systemOperator->SetBlock(0,0, m_systemBlock11);
    m_systemOperator->SetBlock(0,1, m_systemBlock12);
    m_systemOperator->SetBlock(1,0, m_systemBlock21);
    m_systemOperator->SetBlock(1,1, m_systemBlock22);
}

// Builds the linear system matrix from the blocks as a monolithic matrix
void heat::LsqXtFem
:: rebuildSystemMatrix()
{
    resetSystemOperators();

    rebuildSystemBlocks();

    auto systemBlockMatrix
            = new BlockMatrix(m_blockOffsets);
    systemBlockMatrix->SetBlock(0,0, m_systemBlock11);
    systemBlockMatrix->SetBlock(0,1, m_systemBlock12);
    systemBlockMatrix->SetBlock(1,0, m_systemBlock21);
    systemBlockMatrix->SetBlock(1,1, m_systemBlock22);

    m_systemMatrix = systemBlockMatrix->CreateMonolithic();
    if (!m_systemMatrix->ColumnsAreSorted()) {
        m_systemMatrix->SortColumnIndices();
    }
    delete systemBlockMatrix;

    applyBCs(*m_systemMatrix);
}

void heat::LsqXtFem
:: buildSystemMatrix()
{
    buildSystemBlocks();

    auto systemBlockMatrix
            = new BlockMatrix(m_blockOffsets);
    systemBlockMatrix->SetBlock(0,0, m_systemBlock11);
    systemBlockMatrix->SetBlock(0,1, m_systemBlock12);
    systemBlockMatrix->SetBlock(1,0, m_systemBlock21);
    systemBlockMatrix->SetBlock(1,1, m_systemBlock22);

    m_systemMatrix = systemBlockMatrix->CreateMonolithic();
    if (!m_systemMatrix->ColumnsAreSorted()) {
        m_systemMatrix->SortColumnIndices();
    }
    delete systemBlockMatrix;

    applyBCs(*m_systemMatrix);
}

// Builds system matrix from the system operator
// stores only the upper triangle
void heat::LsqXtFem
:: rebuildUpperTriangleOfSystemMatrix()
{
    resetSystemOperators();

    rebuildSystemBlocks();

    auto systemBlockMatrix
            = new BlockMatrix(m_blockOffsets);
    systemBlockMatrix->SetBlock(0,0, m_systemBlock11);
    systemBlockMatrix->SetBlock(0,1, m_systemBlock12);
    systemBlockMatrix->SetBlock(1,0, m_systemBlock21);
    systemBlockMatrix->SetBlock(1,1, m_systemBlock22);

    auto systemMatrixBuf = systemBlockMatrix->CreateMonolithic();
    delete systemBlockMatrix;

    applyBCs(*systemMatrixBuf);

    m_systemMatrix = &getUpperTriangle(*systemMatrixBuf);
    if (!m_systemMatrix->ColumnsAreSorted()) {
        m_systemMatrix->SortColumnIndices();
    }
    delete systemMatrixBuf;
}

void heat::LsqXtFem
:: buildUpperTriangleOfSystemMatrix()
{
    buildSystemBlocks();

    auto systemBlockMatrix
            = new BlockMatrix(m_blockOffsets);
    systemBlockMatrix->SetBlock(0,0, m_systemBlock11);
    systemBlockMatrix->SetBlock(0,1, m_systemBlock12);
    systemBlockMatrix->SetBlock(1,0, m_systemBlock21);
    systemBlockMatrix->SetBlock(1,1, m_systemBlock22);

    auto systemMatrixBuf = systemBlockMatrix->CreateMonolithic();
    delete systemBlockMatrix;

    applyBCs(*systemMatrixBuf);

    m_systemMatrix = &getUpperTriangle(*systemMatrixBuf);
    if (!m_systemMatrix->ColumnsAreSorted()) {
        m_systemMatrix->SortColumnIndices();
    }
    delete systemMatrixBuf;
}

void heat::LsqXtFem
:: rebuildSystemBlocks()
{
    if (m_firstPassOfAssembleSystemBlocks) {
        buildSystemBlock22();
    }

    resetMediumDependentSystemBlocks();
    rebuildSystemBlock11();
    rebuildSystemBlock12AndBlock21();

    if (m_firstPassOfAssembleSystemBlocks) {
        m_firstPassOfAssembleSystemBlocks = false;
    }
}

void heat::LsqXtFem
:: buildSystemBlocks()
{
    buildSystemBlock22();
    buildSystemBlock11();
    buildSystemBlock12AndBlock21();
}

void heat::LsqXtFem
:: buildSystemBlock22()
{
    auto tmp = Add(*m_spatialMass2, *m_spatialStiffness2);
    m_systemBlock22 = OuterProduct(*m_temporalMass, *tmp);
    delete tmp;
}

void heat::LsqXtFem
:: rebuildSystemBlock11()
{
    if (m_firstPassOfAssembleSystemBlocks) {
        buildMediumIndependentSystemBlock11();
    }
    buildMediumDependentSystemBlock11();
}

void heat::LsqXtFem
:: buildSystemBlock11()
{
    buildMediumIndependentSystemBlock11();
    buildMediumDependentSystemBlock11();
}

void heat::LsqXtFem
:: buildMediumIndependentSystemBlock11()
{
    auto tmp = Add(*m_temporalInitial, *m_temporalStiffness);
    m_materialIndependentSystemBlock11
            = OuterProduct(*tmp, *m_spatialMass1);
    delete tmp;
}

void heat::LsqXtFem
:: buildMediumDependentSystemBlock11()
{
    auto tmp = OuterProduct(*m_temporalMass, *m_spatialStiffness1);
    m_systemBlock11 = Add(*m_materialIndependentSystemBlock11, *tmp);
    delete tmp;
}

void heat::LsqXtFem
:: rebuildSystemBlock12AndBlock21()
{
    if (m_firstPassOfAssembleSystemBlocks) {
        buildMediumIndependentSystemBlock12();
    }
    buildMediumDependentSystemBlock12();
    m_systemBlock21 = Transpose(*m_systemBlock12);
}

void heat::LsqXtFem
:: buildSystemBlock12AndBlock21()
{
    buildMediumIndependentSystemBlock12();
    buildMediumDependentSystemBlock12();
    m_systemBlock21 = Transpose(*m_systemBlock12);
}

void heat::LsqXtFem
:: buildMediumIndependentSystemBlock12()
{
    auto tmp = Transpose(*m_temporalGradient);
    m_materialIndependentSystemBlock12
            = OuterProduct(*tmp, *m_spatialDivergence);
    delete tmp;
}

void heat::LsqXtFem
:: buildMediumDependentSystemBlock12()
{
    auto tmp2 = Transpose(*m_spatialGradient);
    auto tmp1 = OuterProduct(*m_temporalMass, *tmp2);
    m_systemBlock12 = Add(-1, *tmp1,
                          -1, *m_materialIndependentSystemBlock12);
    delete tmp1;
    delete tmp2;
}

void heat::LsqXtFem
:: assembleRhs(BlockVector* B) const
{    
    assembleICs(B->GetBlock(0));
    assembleSource(*B);
    applyBCs(*B);
}

// Assembles initial conditions
void heat::LsqXtFem
:: assembleICs(Vector& b) const
{
    int xdimV = m_spatialFeSpaces[0]->GetTrueVSize();
    Vector bx(xdimV);
    assembleSpatialICs(bx);

    Array<int> vdofs;
    const FiniteElement *fe = nullptr;
    Vector elvec, shape;

    int n=0; // only for element (0, dt)
    {
        m_temporalFeSpace->GetElementVDofs (n, vdofs);
        assert(vdofs[0] == 0);

        fe = m_temporalFeSpace->GetFE(n);
        int tNdofs = fe->GetDof();
        shape.SetSize(tNdofs);
        elvec.SetSize(tNdofs*xdimV);

        IntegrationPoint ip;
        ip.Set1w(0.0, 1.0);
        fe->CalcShape(ip, shape);

        for (int j=0; j<tNdofs; j++) {
            for (int k=0; k<xdimV; k++) {
                elvec(k + j*xdimV) = shape(j)*bx(k);
            }
        }
        addVector(vdofs, elvec, b);
    }
}

void heat::LsqXtFem
:: assembleSpatialICs(Vector& b) const
{
    LinearForm ics_form(m_spatialFeSpaces[0]);
    heat::InitialTemperatureCoeff f(m_testCase);
    ics_form.AddDomainIntegrator
            (new DomainLFIntegrator(f));
    ics_form.Assemble();
    b = ics_form.GetData();
}

void heat::LsqXtFem
:: assembleSource(BlockVector& B) const
{
    assembleSourceWithTemporalGradientOfTemperatureBasis(B.GetBlock(0));
    assembleSourceWithSpatialDivergenceOfHeatFluxBasis(B.GetBlock(1));
}

void heat::LsqXtFem
:: assembleSourceWithTemporalGradientOfTemperatureBasis(Vector &b) const
{
    int Nt = m_temporalFeSpace->GetNE();
    int xdimV = m_spatialFeSpaces[0]->GetTrueVSize();
    Vector bx(xdimV);

    Array<int> vdofs;
    ElementTransformation *trans = nullptr;
    const FiniteElement *fe = nullptr;
    Vector elvec;
    DenseMatrix dshape;

    for (int n=0; n<Nt; n++)
    {
        m_temporalFeSpace->GetElementVDofs (n, vdofs);
        trans = m_temporalFeSpace->GetElementTransformation(n);

        // compute elvec
        fe = m_temporalFeSpace->GetFE(n);
        int tNdofs = fe->GetDof();
        dshape.SetSize(tNdofs, 1); // time dimension is 1d
        elvec.SetSize(tNdofs*xdimV);

        int order = 2*fe->GetOrder()+2;
        const IntegrationRule *ir
                = &IntRules.Get(fe->GetGeomType(), order);

        elvec = 0.0;
        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir->IntPoint(i);
            trans->SetIntPoint(&ip);
            fe->CalcPhysDShape(*trans, dshape);

            Vector t;
            trans->Transform(ip, t);
            assembleSpatialSourceWithTemperatureBasis(bx, t(0));

            // elvec += w*shape*bx
            double w = ip.weight*trans->Weight();
            for (int j=0; j<tNdofs; j++){
                double coeff = w*dshape(j,0);
                for (int k=0; k<xdimV; k++) {
                    elvec(k + j*xdimV) += coeff*bx(k);
                }
            }
        }
        addVector(vdofs, elvec, b);
    }
}

void heat::LsqXtFem
:: assembleSpatialSourceWithTemperatureBasis(Vector& b, double t) const
{
    LinearForm source_form(m_spatialFeSpaces[0]);
    heat::SourceCoeff f(m_testCase);
    f.SetTime(t);
    source_form.AddDomainIntegrator
            (new DomainLFIntegrator(f, 2, 4));
    source_form.Assemble();
    b = source_form.GetData();
}

void heat::LsqXtFem
:: assembleSourceWithSpatialDivergenceOfHeatFluxBasis(Vector& b) const
{
    int Nt = m_temporalFeSpace->GetNE();
    int xdimR = m_spatialFeSpaces[1]->GetTrueVSize();
    Vector bx(xdimR);

    Array<int> vdofs;
    ElementTransformation *trans = nullptr;
    const FiniteElement *fe = nullptr;
    Vector elvec, shape;

    for (int n=0; n<Nt; n++)
    {
        m_temporalFeSpace->GetElementVDofs (n, vdofs);
        trans = m_temporalFeSpace->GetElementTransformation(n);

        // compute elvec
        fe = m_temporalFeSpace->GetFE(n);
        int tNdofs = fe->GetDof();
        shape.SetSize(tNdofs);
        elvec.SetSize(tNdofs*xdimR);

        int order = 2*fe->GetOrder()+2;
        const IntegrationRule *ir
                = &IntRules.Get(fe->GetGeomType(), order);

        elvec = 0.0;
        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir->IntPoint(i);
            trans->SetIntPoint(&ip);
            fe->CalcShape(ip, shape);

            Vector t;
            trans->Transform(ip, t);
            assembleSpatialSourceWithSpatialDivergenceOfHeatFluxBasis
                    (bx, t(0));

            // elvec += w*shape*bx
            double w = ip.weight*trans->Weight();
            for (int j=0; j<tNdofs; j++){
                double coeff = w*shape(j);
                for (int k=0; k<xdimR; k++) {
                    elvec(k + j*xdimR) += coeff*bx(k);
                }
            }
        }
        addVector(vdofs, elvec, b);
    }
}

// Adds a local elvec to the global vector b
void heat::LsqXtFem
:: addVector(const Array<int> &vdofs, const Vector& elvec,
             Vector& b) const
{
    int ndofs = elvec.Size();
    int tNdofs = vdofs.Size();
    int xNdofs = ndofs/tNdofs;

    const bool use_dev
            = vdofs.UseDevice() || elvec.UseDevice();
    auto d_elvec = elvec.Read(use_dev);
    auto d_b = b.ReadWrite(use_dev);

    for(int i=0; i<tNdofs; i++)
    {
        const int j = vdofs[i];
        const int ii = i*xNdofs;
        if (j >= 0) {
            const int jj = j*xNdofs;
            for (int k=0; k<xNdofs; k++)
                d_b[k + jj] += d_elvec[k + ii];
        }
        else {
            const int jj = (-1-j)*xNdofs;
            for (int k=0; k<xNdofs; k++)
                d_b[k + jj] -= d_elvec[k + ii];
        }
    }
}

// Applies BCs to SparseMatrix
void heat::LsqXtFem :: applyBCs(SparseMatrix& A) const
{
    for (int k=0; k<m_essentialDofs.Size(); k++) {
        A.EliminateRowCol(m_essentialDofs[k]);
    }
}

// Applies BCs to BlockVector
void heat::LsqXtFem :: applyBCs(BlockVector& B) const
{
    Vector& B0 = B.GetBlock(0);
    B0.SetSubVector(m_essentialDofs, 0.0);
}

// End of file
