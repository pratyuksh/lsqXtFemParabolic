#include "temporal_assembly.hpp"


sparseHeat::TemporalBlockMatrixAssembler
:: TemporalBlockMatrixAssembler(double T, int minLevel, int maxLevel)
    : m_endTime(T),
      m_minLevel(minLevel), m_maxLevel(maxLevel)
{
    m_numLevels = maxLevel - minLevel + 1;
    initialize();
}

void sparseHeat::TemporalBlockMatrixAssembler
:: initialize()
{
    evalMeshCharacteristics();
    evalBlockSizesAndOffsets();
    allocateBlockMatrix();
}

void sparseHeat::TemporalBlockMatrixAssembler
:: evalMeshCharacteristics()
{
    m_numMeshElements.SetSize(m_numLevels);
    m_meshSizes.SetSize(m_numLevels);

    for(int i=0; i<m_numLevels; i++)
    {
        m_numMeshElements[i]
                = static_cast<int>(std::pow(2, m_minLevel+i));
        m_meshSizes[i] = m_endTime / m_numMeshElements[i];
    }
}

void sparseHeat::TemporalBlockMatrixAssembler
:: evalBlockSizesAndOffsets()
{
    m_blockSizes.SetSize(m_numLevels);
    m_blockOffsets.SetSize(m_numLevels+1);

    m_blockOffsets[0] = 0;
    // for the minimum level, the hierarchical basis functions are the
    // same as the standard basis functions
    m_blockSizes[0] = 1+static_cast<int>(std::pow(2, m_minLevel));
    m_blockOffsets[1] = m_blockSizes[0];
    for (int i=1; i<m_numLevels; i++)
    {
        m_blockSizes[i]
                = static_cast<int>(std::pow(2, m_minLevel+i-1));
        m_blockOffsets[i+1] = m_blockSizes[i];
    }
    m_blockOffsets.PartialSum();
}

void sparseHeat::TemporalBlockMatrixAssembler
:: allocateBlockMatrix()
{
    m_blockMatrix = std::make_shared<BlockMatrix>(m_blockOffsets);
    for(int i=0; i<m_numLevels; i++) {
        for(int j=0; j<m_numLevels; j++) {
            SparseMatrix *mat
                    = new SparseMatrix(m_blockSizes[i], m_blockSizes[j]);
            m_blockMatrix->SetBlock(i, j, mat);
        }
    }
    // owns the allocated memory
    m_blockMatrix->owns_blocks = true;
}


void sparseHeat::TemporalCanonicalMatrixAssembler
:: assemble()
{
    auto& mat = m_blockMatrix->GetBlock(0, 0);
    mat.Set(0, 0, 1.0);
    m_blockMatrix->Finalize();
}


void sparseHeat::TemporalMassMatrixAssembler
:: assemble()
{
    // Diagonal blocks
    { // coarsest level will have standard basis functions
        auto& mat = m_blockMatrix->GetBlock(0, 0);
        double h = m_meshSizes[0];

        double hDiv3 = h / 3;
        double hDiv6 = hDiv3 / 2;
        double hTimes2Div3 = 2 * hDiv3;

        for (int i=0; i<m_blockSizes[0]; i++)
        {
            // diagonal element
            if (i == 0 || i == m_blockSizes[0]-1) {
                mat.Set(i, i, hDiv3);
            }
            else {
                mat.Set(i, i, hTimes2Div3);
            }

            // off-diagonal elements
            if (i > 0) {
                mat.Set(i, i-1, hDiv6);
            }
            if (i < m_blockSizes[0]-1) {
                mat.Set(i, i+1, hDiv6);
            }
        }
    }
    // hierarchical basis function after the coarsest level
    for (int m=1; m<m_numLevels; m++)
    {
        auto& mat = m_blockMatrix->GetBlock(m, m);
        double h = m_meshSizes[m];
        double hTimes2Div3 = 2 * h / 3;

        for (int i=0; i<m_blockSizes[m]; i++) {
            mat.Set(i, i, hTimes2Div3);
        }
    }

    // Lower-diagonal blocks
    for (int m=1; m<m_numLevels; m++)
    {
        double hm = m_meshSizes[m];
        for (int n=m-1; n>=0; n--)
        {
            auto& mat = m_blockMatrix->GetBlock(m, n);
            double hn = m_meshSizes[n];

            int factor = static_cast<int>(std::pow(2, m-n));
            int factorDiv2 = factor / 2;
            double hmSquaredDivHn = hm * hm / hn;

            if (n > 0) {
                for (int i=0; i<m_blockSizes[m]; i++)
                {
                    int j = i / factor;
                    int jj = i % factor;
                    jj = {jj >= factorDiv2 ? (factor - jj - 1) : jj};
                    mat.Set(i, j, (2 * jj + 1) * hmSquaredDivHn);
//                    std::cout << m << "\t" << n << "\t" << i << "\t" << jj << std::endl;
                }
            }
            else { // coarsest-level in the column blocks
                for (int i=0; i<m_blockSizes[m]; i++)
                {
                    int jj = i % factor;

                    int jOdd = 2 * (i / factor) + 1;
                    int jjOdd = {jj >= factorDiv2 ? (factor - jj - 1) : jj};

                    int jEven = 2 * (i / factor);
                    jEven = {jj >= factorDiv2 ? jEven+2 : jEven};
                    int jjEven = {jj >= factorDiv2 ? jj - factorDiv2 : factorDiv2 - jj - 1};

                    mat.Set(i, jOdd, (2*jjOdd+1)*hmSquaredDivHn);
                    mat.Set(i, jEven, (2*jjEven+1)*hmSquaredDivHn);
                    // std::cout << m << "\t" << n << "\t" << factor << "\t" << i << "\t" << jOdd << "\t" << jjOdd << std::endl;
                    // std::cout << m << "\t" << n << "\t" << factor << "\t" << i << "\t" << jEven << "\t" << jjEven << std::endl;
                }
            }
        }
    }
    m_blockMatrix->Finalize();

    // Upper-diagonal blocks
    for (int n=1; n<m_numLevels; n++)
        for (int m=0; m<n; m++)
        {
            auto matT = Transpose(m_blockMatrix->GetBlock(n,m));
            m_blockMatrix->GetBlock(m,n)
                    = *matT;
            delete matT;
        }
}


void sparseHeat::TemporalStiffnessMatrixAssembler
:: assemble()
{
    // Diagonal blocks
    { // coarsest level will have standard basis functions
        auto& mat = m_blockMatrix->GetBlock(0, 0);
        double h = m_meshSizes[0];

        double invH = 1./h;
        double invHTimes2 = invH * 2;

        for (int i=0; i<m_blockSizes[0]; i++)
        {
            // diagonal element
            if (i == 0 || i == m_blockSizes[0]-1) {
                mat.Set(i, i, invH);
            }
            else {
                mat.Set(i, i, invHTimes2);
            }

            // off-diagonal elements
            if (i > 0) {
                mat.Set(i, i-1, -invH);
            }
            if (i < m_blockSizes[0]-1) {
                mat.Set(i, i+1, -invH);
            }
        }
    }
    // hierarchical basis function after the coarsest level
    for (int m=1; m<m_numLevels; m++)
    {
        auto& mat = m_blockMatrix->GetBlock(m, m);
        double h = m_meshSizes[m];

        double invH = 1./h;
        double invHTimes2 = invH * 2;

        for (int i=0; i<m_blockSizes[m]; i++) {
            mat.Set(i, i, invHTimes2);
        }
    }
    m_blockMatrix->Finalize();
}


void sparseHeat::TemporalGradientMatrixAssembler
:: assemble()
{
    // Diagonal blocks
    { // coarsest level will have standard basis functions
        auto& mat = m_blockMatrix->GetBlock(0, 0);
        double half = 1./2;

        for (int i=0; i<m_blockSizes[0]; i++)
        {
            // diagonal element
            if (i == 0) {
                mat.Set(i, i, -half);
            }
            else if (i == m_blockSizes[0]-1) {
                mat.Set(i, i, +half);
            }

            // off-diagonal elements
            if (i > 0) {
                mat.Set(i, i-1, -half);
            }
            if (i < m_blockSizes[0]-1) {
                mat.Set(i, i+1, +half);
            }
        }
    }
    // Lower-diagonal blocks
    for (int m=1; m<m_numLevels; m++)
    {
        double hm = m_meshSizes[m];
        for (int n=m-1; n>=0; n--)
        {
            auto& mat = m_blockMatrix->GetBlock(m, n);
            double hn = m_meshSizes[n];

            int factor = static_cast<int>(std::pow(2, m-n));
            int factorDiv2 = factor/2;
            double hmDivHn = hm / hn;

            if (n > 0) {
                for (int i=0; i<m_blockSizes[m]; i++)
                {
                    int j = i / factor;
                    int sign = {i % factor < factorDiv2 ? +1 : -1};
                    mat.Set(i, j, sign*hmDivHn);
                }
            }
            else { // coarsest-level in the column blocks
                for (int i=0; i<m_blockSizes[m]; i++)
                {
                    int jj = i % factor;

                    int jOdd = 2 * (i / factor) + 1;
                    int sign = {jj < factorDiv2 ? +1 : -1};

                    int jEven = 2 * (i / factor);
                    jEven = {jj >= factorDiv2 ? jEven+2 : jEven};

                    mat.Set(i, jOdd, sign*hmDivHn);
                    mat.Set(i, jEven, -sign*hmDivHn);
//                    std::cout << m << "\t" << n << "\t" << factor << "\t" << i << "\t" << jOdd << "\t" << sign << std::endl;
                }
            }
        }
    }
    // Upper-diagonal blocks
    for (int n=1; n<m_numLevels; n++)
    {
        double hn = m_meshSizes[n];
        for (int m=n-1; m>=0; m--)
        {
            auto& mat = m_blockMatrix->GetBlock(m, n);
            double hm = m_meshSizes[m];

            int factor = static_cast<int>(std::pow(2, n-m));
            int factorDiv2 = factor/2;
            double hnDivHm = hn / hm;

            if (m > 0) {
                for (int j=0; j<m_blockSizes[n]; j++)
                {
                    int i = j / factor;
                    int sign = {j % factor < factorDiv2 ? -1 : +1};
                    mat.Set(i, j, sign*hnDivHm);
                }
            }
            else { // coarsest-level in the column blocks
                for (int j=0; j<m_blockSizes[n]; j++)
                {
                    int ii = j % factor;

                    int iOdd = 2 * (j / factor) + 1;
                    int sign = {ii < factorDiv2 ? -1 : +1};

                    int iEven = 2 * (j / factor);
                    iEven = {ii >= factorDiv2 ? iEven+2 : iEven};

                    mat.Set(iOdd, j, sign*hnDivHm);
                    mat.Set(iEven, j, -sign*hnDivHm);
                }
            }
        }
    }
    m_blockMatrix->Finalize();
}


