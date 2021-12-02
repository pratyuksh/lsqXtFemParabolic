#include "coefficients.hpp"

using namespace mfem;


//-----------------------//
//  Medium Coefficients  //
//----------------------//

// Heat medium scalar coefficient
double heat::MediumScalarCoeff
:: Eval (ElementTransformation& T,
         const IntegrationPoint& ip)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    return m_testCase->medium(transip);
}

// Heat medium tensor coefficient
void heat::MediumTensorCoeff
:: Eval (DenseMatrix& M, ElementTransformation& T,
         const IntegrationPoint& ip)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    M = m_testCase->mediumTensor(transip);
}


//-------------------------------//
//  Exact Solution Coefficients  //
//-------------------------------//

// Heat exact temperature coefficient
double heat::ExactTemperatureCoeff
:: Eval (ElementTransformation& T,
         const IntegrationPoint& ip)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    return (m_testCase->temperatureSol(transip, GetTime()));
}

double heat::ExactTemperatureCoeff
:: Eval (ElementTransformation& T,
         const IntegrationPoint& ip,
         double t)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    return (m_testCase->temperatureSol(transip, t));
}

// Heat exact flux coefficient
void heat::ExactHeatFluxCoeff
:: Eval (Vector& v,
         ElementTransformation& T,
         const IntegrationPoint& ip)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    v = m_testCase->heatFluxSol(transip, GetTime());
}

void heat::ExactHeatFluxCoeff
:: Eval (Vector& v,
         ElementTransformation& T,
         const IntegrationPoint& ip,
         double t)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    v = m_testCase->heatFluxSol(transip, t);
}

// Heat exact temperature spatial gradient coefficient
void heat::ExactTemperatureSpatialGradCoeff
:: Eval (Vector& v,
         ElementTransformation& T,
         const IntegrationPoint& ip)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    v = m_testCase->temperatureSpatialGradientSol(transip, GetTime());
}

void heat::ExactTemperatureSpatialGradCoeff
:: Eval (Vector& v,
         ElementTransformation& T,
         const IntegrationPoint& ip,
         double t)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    v = m_testCase->temperatureSpatialGradientSol(transip, t);
}

// Heat exact temperature time-gradient coefficient
double heat::ExactTemperatureTemporalGradCoeff
:: Eval (ElementTransformation& T,
         const IntegrationPoint& ip)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    return (m_testCase->temperatureTemporalGradientSol
            (transip, GetTime()));
}

double heat::ExactTemperatureTemporalGradCoeff
:: Eval (ElementTransformation& T,
         const IntegrationPoint& ip,
         double t)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    return (m_testCase->temperatureTemporalGradientSol
            (transip, t));
}

// Heat Laplacian coefficient
double heat::LaplacianCoeff
:: Eval (ElementTransformation& T,
         const IntegrationPoint& ip)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    return (m_testCase->laplacian(transip, GetTime()));
}

double heat::LaplacianCoeff
:: Eval (ElementTransformation& T,
         const IntegrationPoint& ip,
         double t)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    return (m_testCase->laplacian(transip, t));
}

//---------------------------------//
//  Initial Solution Coefficients  //
//---------------------------------//

// Heat initial temperature coefficient
double heat::InitialTemperatureCoeff
:: Eval (ElementTransformation& T,
         const IntegrationPoint& ip)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    return (m_testCase->initTemperature(transip));
}


//------------------------------------//
//  Boundary Temperature Coefficient  //
//------------------------------------//

// Heat boundary temperature coefficient
double heat::BdryTemperatureCoeff
:: Eval (ElementTransformation& T,
         const IntegrationPoint& ip)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    return (m_testCase->bdryTemperature(transip, GetTime()));
}

double heat::BdryTemperatureCoeff
:: Eval (ElementTransformation& T,
         const IntegrationPoint& ip,
         double t)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    return (m_testCase->bdryTemperature(transip, t));
}

//-----------------------//
//  Source Coefficients  //
//-----------------------//

// Heat source coefficient
double heat::SourceCoeff
:: Eval (ElementTransformation& T,
         const IntegrationPoint& ip)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    return m_testCase->source(transip, GetTime());
}

double heat::SourceCoeff
:: Eval (ElementTransformation& T,
         const IntegrationPoint& ip,
         double t)
{
    double x[3];
    Vector transip(x, 3);
    T.Transform(ip, transip);
    return m_testCase->source(transip, t);
}


// End of file
