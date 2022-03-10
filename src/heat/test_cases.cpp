#include "test_cases.hpp"

using namespace mfem;

// Dummy
// Value 1 everywhere, except boundary
// Homogeneous Dirichlet BCs
// Zero source
double heat::TestCase <Dummy>
:: temperatureSol(const Vector& x, const double t) const
{
    return 1;
}

Vector heat::TestCase <Dummy>
:: heatFluxSol (const Vector& x, const double t) const
{
    Vector q(m_dim);
    q = 0.;
    return q;
}

Vector heat::TestCase <Dummy>
:: temperatureSpatialGradientSol (const Vector& x, const double t) const
{
    Vector gradu(m_dim);
    gradu = 0.;
    return gradu;
}

double heat::TestCase <Dummy>
:: temperatureTemporalGradientSol
(const Vector& x, const double t) const
{
    return 0;
}

double heat::TestCase <Dummy>
:: source (const Vector&, const double) const
{
    return 0;
}

void heat::TestCase <Dummy>
:: setBdryDirichlet(Array<int>& bdr_marker) const
{
    bdr_marker = 1;
}

// UnitSquareTest1, Smooth solution
// Non-zero ICs
// Homogeneous Dirichlet BCs
// Zero source
double heat::TestCase <UnitSquareTest1>
:: medium(const Vector&) const
{
    return std::exp(m_rvar);
}

DenseMatrix heat::TestCase <UnitSquareTest1>
:: mediumTensor(const Vector& x) const
{
    DenseMatrix med(m_dim);
    med(0,0) = med(1,1) = 1;
    med(0,1) = med(1,0) = 0;
    med *= medium(x);
    return med;
}

double heat::TestCase <UnitSquareTest1>
:: temperatureSol(const Vector& x, const double t) const
{
    double coeff = medium(x);
    double temperature
            = sin(M_PI*x(0))*sin(M_PI*x(1))
            *exp(-2*M_PI*M_PI*coeff*t);
    return temperature;
}

Vector heat::TestCase <UnitSquareTest1>
:: heatFluxSol (const Vector& x, const double t) const
{
    return this->temperatureSpatialGradientSol(x, t);
}

Vector heat::TestCase <UnitSquareTest1>
:: temperatureSpatialGradientSol (const Vector& x, const double t) const
{
    double coeff = medium(x);
    Vector dudx(m_dim);
    dudx(0) = M_PI*cos(M_PI*x(0))*sin(M_PI*x(1));
    dudx(1) = M_PI*sin(M_PI*x(0))*cos(M_PI*x(1));
    dudx *= exp(-2*M_PI*M_PI*coeff*t);
    return dudx;
}

double heat::TestCase <UnitSquareTest1>
:: temperatureTemporalGradientSol
(const Vector& x, const double t) const
{
    double coeff = medium(x);
    double dudt = sin(M_PI*x(0))*sin(M_PI*x(1))
            *(-2*M_PI*M_PI*coeff)
            *exp(-2*M_PI*M_PI*coeff*t);
    return dudt;
}

double heat::TestCase <UnitSquareTest1>
:: source (const Vector&, const double) const
{
    double f = 0;
    return f;
}

void heat::TestCase <UnitSquareTest1>
:: setBdryDirichlet(Array<int>& bdr_marker) const
{
    bdr_marker = 1;
}

// UnitSquareTest2, Smooth solution
// Non-zero ICs
// Homogeneous Dirichlet BCs
// Non-zero source
double heat::TestCase <UnitSquareTest2>
:: temperatureSol(const Vector& x, const double t) const
{
    double temperature
            = sin(M_PI*x(0))*sin(M_PI*x(1))*cos(M_PI*t);
    return temperature;
}

Vector heat::TestCase <UnitSquareTest2>
:: heatFluxSol (const Vector& x, const double t) const
{
    return temperatureSpatialGradientSol(x, t);
}

Vector heat::TestCase <UnitSquareTest2>
:: temperatureSpatialGradientSol (const Vector& x, const double t) const
{
    Vector dudx(m_dim);
    dudx(0) = M_PI*cos(M_PI*x(0))*sin(M_PI*x(1));
    dudx(1) = M_PI*sin(M_PI*x(0))*cos(M_PI*x(1));
    dudx *= cos(M_PI*t);
    return dudx;
}

double heat::TestCase <UnitSquareTest2>
:: temperatureTemporalGradientSol
(const Vector& x, const double t) const
{
    double dudt = sin(M_PI*x(0))*sin(M_PI*x(1))
            *(-M_PI)*sin(M_PI*t);
    return dudt;
}

double heat::TestCase <UnitSquareTest2>
:: source (const Vector& x, const double t) const
{
    double f = M_PI*sin(M_PI*x(0))*sin(M_PI*x(1));
    f *= (2*M_PI*cos(M_PI*t) - sin(M_PI*t));
    return f;
}

void heat::TestCase <UnitSquareTest2>
:: setBdryDirichlet(Array<int>& bdr_marker) const
{
    bdr_marker = 1;
}

// UnitSquareTest3, Smooth solution
// Zero ICs
// Homogeneous Dirichlet BCs
// Non-zero source
double heat::TestCase <UnitSquareTest3>
:: medium(const Vector& x) const
{
    return std::move(perturb(x));
}

DenseMatrix heat::TestCase <UnitSquareTest3>
:: mediumTensor(const Vector& x) const
{
    DenseMatrix med(m_dim);
    med(0,0) = med(1,1) = 1;
    med(0,1) = med(1,0) = 0;
    med *= medium(x);
    return med;
}

double heat::TestCase <UnitSquareTest3>
:: temperatureSol(const Vector& x, const double t) const
{
    double temperature
            = sin(M_PI*x(0))*sin(M_PI*x(1))*sin(M_PI*t);
    return temperature;
}

Vector heat::TestCase <UnitSquareTest3>
:: heatFluxSol (const Vector& x, const double t) const
{
    return temperatureSpatialGradientSol(x, t);
}

Vector heat::TestCase <UnitSquareTest3>
:: temperatureSpatialGradientSol (const Vector& x, const double t) const
{
    Vector dudx(m_dim);
    dudx(0) = M_PI*cos(M_PI*x(0))*sin(M_PI*x(1));
    dudx(1) = M_PI*sin(M_PI*x(0))*cos(M_PI*x(1));
    dudx *= sin(M_PI*t);
    return dudx;
}

double heat::TestCase <UnitSquareTest3>
:: temperatureTemporalGradientSol
(const Vector& x, const double t) const
{
    double dudt = sin(M_PI*x(0))*sin(M_PI*x(1))
            *M_PI*cos(M_PI*t);
    return dudt;
}

double heat::TestCase <UnitSquareTest3>
:: source (const Vector& x, const double t) const
{
    double f = M_PI*sin(M_PI*x(0))*sin(M_PI*x(1));
    f *= (2*M_PI*sin(M_PI*t) + cos(M_PI*t));
    return f;
}

void heat::TestCase <UnitSquareTest3>
:: setBdryDirichlet(Array<int>& bdr_marker) const
{
    bdr_marker = 1;
}

double heat::TestCase <UnitSquareTest3>
:: perturb (const Vector& x) const
{
    double val = std::exp(m_rvar*(sin(M_PI*x(0))
                                  + sin(M_PI*x(1))));
    return val;
}

// UnitSquareTest4, Smooth solution
// Zero ICs
// Homogeneous Dirichlet BCs
// Zero source
DenseMatrix heat::TestCase <UnitSquareTest4>
:: mediumTensor(const Vector&) const
{
    DenseMatrix med(m_dim);
    med(0,0) = 1./2.;
    med(1,1) = 2./3.;
    med(0,1) = med(1,0) = 1./4.;
    return med;
}

double heat::TestCase <UnitSquareTest4>
:: temperatureSol(const Vector& x, const double t) const
{
    double temperature
            = sin(M_PI*x(0))*sin(M_PI*x(1))*sin(M_PI*t);
    return temperature;
}

Vector heat::TestCase <UnitSquareTest4>
:: heatFluxSol (const Vector& x, const double t) const
{
    Vector q(m_dim);
    auto tmp = temperatureSpatialGradientSol(x, t);
    auto medMat = mediumTensor(x);
    medMat.Mult(tmp, q);
    return q;
}

Vector heat::TestCase <UnitSquareTest4>
:: temperatureSpatialGradientSol (const Vector& x, const double t) const
{
    Vector dudx(m_dim);
    dudx(0) = M_PI*cos(M_PI*x(0))*sin(M_PI*x(1));
    dudx(1) = M_PI*sin(M_PI*x(0))*cos(M_PI*x(1));
    dudx *= sin(M_PI*t);
    return dudx;
}

double heat::TestCase <UnitSquareTest4>
:: temperatureTemporalGradientSol
(const Vector& x, const double t) const
{
    double dudt = sin(M_PI*x(0))*sin(M_PI*x(1))
            *M_PI*cos(M_PI*t);
    return dudt;
}

double heat::TestCase <UnitSquareTest4>
:: source (const Vector& x, const double t) const
{
    double f1 = sin(M_PI*x(0))*sin(M_PI*x(1));
    f1 *= M_PI*cos(M_PI*t);
    double f2 = -7./6.*sin(M_PI*x(0))*sin(M_PI*x(1));
    f2 += 1./2.*cos(M_PI*x(0))*cos(M_PI*x(1));
    f2 *= -M_PI*M_PI*sin(M_PI*t);
    return f1+f2;
}

void heat::TestCase <UnitSquareTest4>
:: setBdryDirichlet(Array<int>& bdr_marker) const
{
    bdr_marker = 1;
}

// PeriodicUnitSquareTest1, Smooth solution
double heat::TestCase <PeriodicUnitSquareTest1>
:: medium(const Vector& x) const
{
    return std::move(perturb(x));
}

DenseMatrix heat::TestCase <PeriodicUnitSquareTest1>
:: mediumTensor(const Vector& x) const
{
    DenseMatrix med(m_dim);
    med(0,0) = med(1,1) = 1;
    med(0,1) = med(1,0) = 0;
    med *= medium(x);
    return med;
}

double heat::TestCase <PeriodicUnitSquareTest1>
:: temperatureSol(const Vector& x, const double t) const
{
    double alpha = 4*M_PI*M_PI;
    double ux = 200*(sin(2*M_PI*x(0)) + sin(2*M_PI*x(1)));
    double ut = (1 - exp(-alpha*t))/alpha;
    double temperature = ux*ut;
    return temperature;
}

Vector heat::TestCase <PeriodicUnitSquareTest1>
:: heatFluxSol (const Vector& x, const double t) const
{
    Vector q(m_dim);
    auto tmp = temperatureSpatialGradientSol(x, t);
    auto medMat = mediumTensor(x);
    medMat.Mult(tmp, q);
    return q;
}

Vector heat::TestCase <PeriodicUnitSquareTest1>
:: temperatureSpatialGradientSol (const Vector& x, const double t) const
{
    Vector dudx(m_dim);
    double alpha = 4*M_PI*M_PI;
    double ut = (1 - exp(-alpha*t))/alpha;
    dudx(0) = 400*M_PI*cos(2*M_PI*x(0));
    dudx(1) = 400*M_PI*cos(2*M_PI*x(1));
    dudx *= ut;
    return dudx;
}

double heat::TestCase <PeriodicUnitSquareTest1>
:: temperatureTemporalGradientSol
(const Vector& x, const double t) const
{
    double alpha = 4*M_PI*M_PI;
    double ux = 200*(sin(2*M_PI*x(0)) + sin(2*M_PI*x(1)));
    double dudt = ux*exp(-alpha*t);
    return dudt;
}

double heat::TestCase <PeriodicUnitSquareTest1>
:: source (const Vector& x, const double) const
{
    double f = 200*(sin(2*M_PI*x(0)) + sin(2*M_PI*x(1)));
    return f;
}

double heat::TestCase <PeriodicUnitSquareTest1>
:: perturb (const Vector& x) const
{
    double val = std::exp(m_rvar*(sin(2*M_PI*x(0))
                                  + sin(2*M_PI*x(1))));
    return val;
}

// LShapedTest1, Singular solution
// Non-zero ICs
// Homogeneous Dirichlet BCs
// Non-zero source
double heat::TestCase <LShapedTest1>
:: radius(const Vector& x) const
{
    double r = x.DistanceTo(m_singularCorner);
    return r;
}

double heat::TestCase <LShapedTest1>
:: polarAngle(const Vector& x) const
{
    double theta = std::atan2(-x[0],x[1])
            + (x[0] >= 0 && x[1] <= 0 ? 2*M_PI : 0);
    return theta;
}

double heat::TestCase <LShapedTest1>
:: temperatureSol(const Vector& x, const double t) const
{
    double r = radius(x);
    double theta = polarAngle(x);
    double temperature
            = pow(r, m_gamma)*sin(m_gamma*theta)
            *(1 - x[0]*x[0])*(1 - x[1]*x[1])
            *cos(M_PI*t);
    return temperature;
}

Vector heat::TestCase <LShapedTest1>
:: heatFluxSol (const Vector& x, const double t) const
{
    return temperatureSpatialGradientSol(x, t);
}

Vector heat::TestCase <LShapedTest1>
:: temperatureSpatialGradientSol (const Vector& x, const double t) const
{
    double r = radius(x);
    double theta = polarAngle(x);

    Vector dudx(m_dim);
    double coeff = m_gamma-1;

    dudx(0) = -m_gamma*pow(r, coeff)*cos(coeff*theta)
            *(1 - x[0]*x[0]);
    dudx(0) += pow(r, m_gamma)*sin(m_gamma*theta)*(- 2*x[0]);
    dudx(0) *= (1. - x[1]*x[1]);

    dudx(1) = m_gamma*pow(r, coeff)*sin(coeff*theta)
            *(1 - x[1]*x[1]);
    dudx(1) += pow(r, m_gamma)*sin(m_gamma*theta)*(- 2*x[1]);
    dudx(1) *= (1 - x[0]*x[0]);

    dudx *= cos(M_PI*t);

    return dudx;
}

double heat::TestCase <LShapedTest1>
:: temperatureTemporalGradientSol
(const Vector& x, const double t) const
{
    double r = radius(x);
    double theta = polarAngle(x);
    double dudt
            = pow(r, m_gamma)*sin(m_gamma*theta)
            *(1./4. - x[0]*x[0])*(1./4. - x[1]*x[1])
            *(-M_PI)*sin(M_PI*t);
    return dudt;
}

double heat::TestCase <LShapedTest1>
:: source (const Vector& x, const double t) const
{
    double r = radius(x);
    double theta = polarAngle(x);
    double coeff = m_gamma-1;

    double f = -M_PI*pow(r, m_gamma)*sin(m_gamma*theta)
            *(1 - x[0]*x[0])*(1 - x[1]*x[1])
            *sin(M_PI*t);

    double f1 = 4*m_gamma*pow(r, coeff)*cos(coeff*theta)*x[0];
    f1 -= 2*pow(r,m_gamma)*sin(m_gamma*theta);
    f1 *= (1 - x[1]*x[1]);

    double f2 = -4*m_gamma*pow(r, coeff)*sin(coeff*theta)*x[1];
    f2 -= 2*pow(r,m_gamma)*sin(m_gamma*theta);
    f2 *= (1 - x[0]*x[0]);

    f -= (f1 + f2)*cos(M_PI*t);

    return f;
}

void heat::TestCase <LShapedTest1>
:: setBdryDirichlet(Array<int>& bdr_marker) const
{
    bdr_marker = 1;
}

// LShapedTest2, Singular solution
// Zero ICs
// Homogeneous Dirichlet BCs
// Non-zero source
double heat::TestCase <LShapedTest2>
:: radius(const Vector& x) const
{
    double r = x.DistanceTo(m_singularCorner);
    return r;
}

double heat::TestCase <LShapedTest2>
:: polarAngle(const Vector& x) const
{
    double theta = std::atan2(-x[0],x[1])
            + (x[0] >= 0 && x[1] <= 0 ? 2*M_PI : 0);
    return theta;
}

double heat::TestCase <LShapedTest2>
:: temperatureSol(const Vector& x, const double t) const
{
    double r = radius(x);
    double theta = polarAngle(x);
    double temperature
            = pow(r, m_gamma)*sin(m_gamma*theta)
            *(1 - x[0]*x[0])*(1 - x[1]*x[1])
            *sin(M_PI*t);
    return temperature;
}

Vector heat::TestCase <LShapedTest2>
:: heatFluxSol (const Vector& x, const double t) const
{
    return temperatureSpatialGradientSol(x, t);
}

Vector heat::TestCase <LShapedTest2>
:: temperatureSpatialGradientSol (const Vector& x, const double t) const
{
    double r = radius(x);
    double theta = polarAngle(x);

    Vector dudx(m_dim);
    double coeff = m_gamma-1;

    dudx(0) = -m_gamma*pow(r, coeff)*cos(coeff*theta)
            *(1 - x[0]*x[0]);
    dudx(0) += pow(r, m_gamma)*sin(m_gamma*theta)*(- 2*x[0]);
    dudx(0) *= (1 - x[1]*x[1]);

    dudx(1) = m_gamma*pow(r, coeff)*sin(coeff*theta)
            *(1 - x[1]*x[1]);
    dudx(1) += pow(r, m_gamma)*sin(m_gamma*theta)*(- 2*x[1]);
    dudx(1) *= (1 - x[0]*x[0]);

    dudx *= sin(M_PI*t);

    return dudx;
}

double heat::TestCase <LShapedTest2>
:: temperatureTemporalGradientSol
(const Vector& x, const double t) const
{
    double r = radius(x);
    double theta = polarAngle(x);
    double temperature
            = pow(r, m_gamma)*sin(m_gamma*theta)
            *(1 - x[0]*x[0])*(1 - x[1]*x[1])
            *M_PI*cos(M_PI*t);
    return temperature;
}

double heat::TestCase <LShapedTest2>
:: source (const Vector& x, const double t) const
{
    double r = radius(x);
    double theta = polarAngle(x);
    double coeff = m_gamma-1;

    double f = M_PI*pow(r, m_gamma)*sin(m_gamma*theta)
            *(1 - x[0]*x[0])*(1 - x[1]*x[1])
            *cos(M_PI*t);

    double f1 = 4*m_gamma*pow(r, coeff)*cos(coeff*theta)*x[0];
    f1 -= 2*pow(r,m_gamma)*sin(m_gamma*theta);
    f1 *= (1 - x[1]*x[1]);

    double f2 = -4*m_gamma*pow(r, coeff)*sin(coeff*theta)*x[1];
    f2 -= 2*pow(r,m_gamma)*sin(m_gamma*theta);
    f2 *= (1 - x[0]*x[0]);

    f -= (f1 + f2)*sin(M_PI*t);

    return f;
}

void heat::TestCase <LShapedTest2>
:: setBdryDirichlet(Array<int>& bdr_marker) const
{
    bdr_marker = 1;
}

// LShapedTest3, Singular solution
// Zero ICs
// Homogeneous Dirichlet BCs
// Constant source 1
double heat::TestCase <LShapedTest3>
:: temperatureSol(const Vector& x, const double t) const {
    return 0;
}

Vector heat::TestCase <LShapedTest3>
:: heatFluxSol (const Vector& x, const double t) const {
    return temperatureSpatialGradientSol(x, t);
}

Vector heat::TestCase <LShapedTest3>
:: temperatureSpatialGradientSol (const Vector& x, const double t) const
{
    Vector dudx(m_dim);
    dudx = 0.;
    return dudx;
}

double heat::TestCase <LShapedTest3>
:: temperatureTemporalGradientSol
(const Vector& x, const double t) const {
    return 0;
}

double heat::TestCase <LShapedTest3>
:: source (const Vector& x, const double t) const {
    return 1;
}

void heat::TestCase <LShapedTest3>
:: setBdryDirichlet(Array<int>& bdr_marker) const {
    bdr_marker = 1;
}

// UnitCubeTest1, Smooth solution
// Non-zero ICs
// Homogeneous Dirichlet BCs
// Non-zero source
double heat::TestCase <UnitCubeTest1>
:: temperatureSol(const Vector& x, const double t) const
{
    double temperature
            = sin(M_PI*x(0))*sin(M_PI*x(1))*sin(M_PI*x(2))*cos(M_PI*t);
    return temperature;
}

Vector heat::TestCase <UnitCubeTest1>
:: heatFluxSol (const Vector& x, const double t) const
{
    return temperatureSpatialGradientSol(x, t);
}

Vector heat::TestCase <UnitCubeTest1>
:: temperatureSpatialGradientSol (const Vector& x, const double t) const
{
    Vector dudx(m_dim);
    dudx(0) = M_PI*cos(M_PI*x(0))*sin(M_PI*x(1))*sin(M_PI*x(2));
    dudx(1) = M_PI*sin(M_PI*x(0))*cos(M_PI*x(1))*sin(M_PI*x(2));
    dudx(2) = M_PI*sin(M_PI*x(0))*sin(M_PI*x(1))*cos(M_PI*x(2));
    dudx *= cos(M_PI*t);
    return dudx;
}

double heat::TestCase <UnitCubeTest1>
:: temperatureTemporalGradientSol
(const Vector& x, const double t) const
{
    double dudt = sin(M_PI*x(0))*sin(M_PI*x(1))*sin(M_PI*x(2))
            *(-M_PI)*sin(M_PI*t);
    return dudt;
}

double heat::TestCase <UnitCubeTest1>
:: source (const Vector& x, const double t) const
{
    double f = M_PI*sin(M_PI*x(0))*sin(M_PI*x(1))*sin(M_PI*x(2));
    f *= (3*M_PI*cos(M_PI*t) - sin(M_PI*t));
    return f;
}

void heat::TestCase <UnitCubeTest1>
:: setBdryDirichlet(Array<int>& bdr_marker) const
{
    bdr_marker = 1;
}

// FicherCubeTest1, Singular solution
// Zero ICs
// Homogeneous Dirichlet BCs
// Constant source 1
double heat::TestCase <FicheraCubeTest1>
:: temperatureSol(const Vector& x, const double t) const {
    return 0;
}

Vector heat::TestCase <FicheraCubeTest1>
:: heatFluxSol (const Vector& x, const double t) const {
    return temperatureSpatialGradientSol(x, t);
}

Vector heat::TestCase <FicheraCubeTest1>
:: temperatureSpatialGradientSol (const Vector& x, const double t) const
{
    Vector dudx(m_dim);
    dudx = 0.;
    return dudx;
}

double heat::TestCase <FicheraCubeTest1>
:: temperatureTemporalGradientSol
(const Vector& x, const double t) const {
    return 0;
}

double heat::TestCase <FicheraCubeTest1>
:: source (const Vector& x, const double t) const {
    return 1;
}

void heat::TestCase <FicheraCubeTest1>
:: setBdryDirichlet(Array<int>& bdr_marker) const {
    bdr_marker = 1;
}

// End of file
