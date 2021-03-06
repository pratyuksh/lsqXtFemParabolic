#include "utilities.hpp"
#include <assert.h>

using namespace mfem;


// Deletes allocated memory and sets pointer to nullptr
template<typename T>
void reset(T* &obj) {
    delete obj;
    obj = nullptr;
}

// free matrix memory and set pointer to null
void clear (SparseMatrix* &mat) {
    delete mat;
    mat = nullptr;
}

// zero function
void zeroFn(const Vector&, Vector& f) {
   f = 0.0;
}

// Returns the upper-triangular, including the diagonal,
// for the input matrix A
SparseMatrix& getUpperTriangle(const SparseMatrix& A)
{
    int numRows = A.NumRows();
    int numCols = A.NumCols();
    const int *iA = A.GetI();
    const int *jA = A.GetJ();
    const double *dA = A.GetData();

    SparseMatrix* UA = new SparseMatrix(numRows, numCols);
    for (int i=0; i<numRows; i++) {
        for(int k=iA[i]; k<iA[i+1]; k++) {
            int j = jA[k];
            if (j >= i) { UA->_Add_(i, j, dA[k]); }
        }
    }
    UA->Finalize();

    return *UA;
}

// Returns cell center of element i of given mesh
void getElementCenter(const std::shared_ptr<Mesh>& mesh,
                      int i, Vector& center)
{
    int geom = mesh->GetElementBaseGeometry(i);
    ElementTransformation *trans
            = mesh->GetElementTransformation(i);
    trans->Transform(Geometries.GetCenter(geom), center);
}



using namespace mymfem;

// Locates a given physical point x to a mesh element
// using the initial guess initElId
std::pair <Array<int>, Array<IntegrationPoint>>
PointLocator :: operator() (const DenseMatrix& X,
                            int& initElId) const
{
    int numPoints = X.NumCols();
    Array<int> elIds(numPoints);
    Array<IntegrationPoint> ips(numPoints);

    Vector x(X.NumRows());
    for (int i=0; i<numPoints; i++)
    {
        X.GetColumn(i, x);
        std::tie (elIds[i], ips[i]) = (*this)(x, initElId);
        initElId = elIds[i];
    }

    return {elIds, ips};
}


// Locates a given physical point x to a mesh element
// using the initial guess initElId
std::pair <int, IntegrationPoint>
PointLocator :: operator() (const Vector& x,
                            const int initElId) const
{
    auto [info_, ip_] = computeRefPoint(x, initElId);
    if (info_ == InverseElementTransformation::Inside) {
        return {initElId, ip_};
    }

    bool found = false;
    int info;
    int elId = initElId, oldElId = initElId;
    IntegrationPoint ip;

    Array <int> vertices;
    Vector z (m_mesh->Dimension());
    double minDist = std::numeric_limits<double>::max();

    while(!found)
    {
        // find the element closest to point x amongst the
        // neighbours of the vertices of the current element
        m_mesh->GetElementVertices(elId, vertices);
        for (int i=0; i < vertices.Size(); i++)
        {
            int v = vertices[i];
            int ne = m_vToEl->RowSize(v);
            const int* els = m_vToEl->GetRow(v);
            //std::cout << "\n\nFor vertex: " << i << std::endl;
            for (int j=0; j<ne; j++)
            {
                if (els[j] == elId) {continue;}

                m_mesh->GetElementTransformation(els[j])
                        ->Transform(Geometries.GetCenter
                         (m_mesh->GetElementBaseGeometry
                          (els[j])), z);
                double dist = z.DistanceTo(x.GetData());

                //std::cout << els[j] << "\t" << elId << "\t"
                //          << dist << "\t"
                //          << minDist << std::endl;
                if (dist < minDist)
                {
                    minDist = dist;
                    elId = els[j];
                    //std::cout << "Minimum: "
                    //          << minDist << "\t"
                    //          << elId << std::endl;
                }
            }
        }

        if (elId == oldElId) {
            // if elId is not upadted, search
            // all neighbours of its vertices
            for (int i=0; i < vertices.Size(); i++)
            {
                int v = vertices[i];
                int ne = m_vToEl->RowSize(v);
                const int* els = m_vToEl->GetRow(v);

                for (int j=0; j<ne; j++)
                {
                    if (els[j] == elId) {continue;}
                    //std::cout << "Search neighbours: "
                    //          << x(0) << "\t"
                    //          << x(1) << "\t"
                    //          << els[j] << std::endl;
                    std::tie (info, ip) = computeRefPoint(x, els[j]);
                    if (info == InverseElementTransformation::Inside) {
                        //std::cout << "Exception: "
                        //          << els[j] << std::endl;
                        return {els[j],ip};
                    }
                }
            }
            // in case the neighbour search fails,
            // loop over other elements
            //if (elId < m_mesh->GetNE()-1) { elId++;}
            //else { elId = 0; }
        }
        else {
            // if elId is upadted, check if it contains point x
            std::tie (info, ip) = computeRefPoint(x, elId);
            if (info == InverseElementTransformation::Inside) {
                found = true;
            }
            oldElId = elId;
        }
    }

    return {elId, ip};
}

// Generates a table of vertices,
// which are shared by different processors.
// Needed when a ParMesh is written to one file,
// the same vertices at the partition interfaces
// will have a different numbering for each processor
void PointLocator :: getSharedVerticesTable() const
{
    int nv = m_mesh->GetNV();
    int dim = m_mesh->Dimension();

    m_sharedVertices->MakeI (nv);
    for (int i=0; i<nv; i++)
    {
        Vector v1(m_mesh->GetVertex(i), dim);
        for (int j=i+1; j<nv; j++)
        {
            if (v1.DistanceTo(m_mesh->GetVertex(j)) < m_TOL) {
                m_sharedVertices->AddAColumnInRow(i);
                m_sharedVertices->AddAColumnInRow(j);
            }
        }
    }
    m_sharedVertices->MakeJ();

    for (int i=0; i<nv; i++)
    {
        Vector v1(m_mesh->GetVertex(i), dim);
        for (int j=i+1; j<nv; j++)
        {
            if (v1.DistanceTo(m_mesh->GetVertex(j)) < m_TOL) {
                m_sharedVertices->AddConnection(i, j);
                m_sharedVertices->AddConnection(j, i);
            }
        }
    }
    m_sharedVertices->ShiftUpI();
}

// Generates the vertex to element table,
// uses the shared vertices table
void PointLocator :: getVerticesToElementsTable() const
{
    int ne = m_mesh->GetNE();
    int nv = m_mesh->GetNV();

    m_vToEl->MakeI(nv);

    Array<int> v, sv;
    for (int i=0; i<ne; i++)
    {
        m_mesh->GetElementVertices(i, v);
        for (int j=0; j<v.Size(); j++) {
            m_vToEl->AddAColumnInRow(v[j]);

            int nsv = m_sharedVertices->RowSize(v[j]);
            m_sharedVertices->GetRow(v[j], sv);
            for (int k=0; k<nsv; k++) {
                m_vToEl->AddAColumnInRow(sv[k]);
            }
        }
    }
    m_vToEl->MakeJ();

    for (int i=0; i<ne; i++)
    {
        m_mesh->GetElementVertices(i, v);
        for (int j=0; j<v.Size(); j++) {
            m_vToEl->AddConnection(v[j], i);

            int nsv = m_sharedVertices->RowSize(v[j]);
            m_sharedVertices->GetRow(v[j], sv);
            for (int k=0; k<nsv; k++) {
                m_vToEl->AddConnection(sv[k], i);
            }
        }
    }
    m_vToEl->ShiftUpI();
}

// End of file
