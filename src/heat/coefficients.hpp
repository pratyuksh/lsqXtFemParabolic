#ifndef HEAT_COEFFICIENTS_HPP
#define HEAT_COEFFICIENTS_HPP

#include "test_cases.hpp"


namespace heat {

//-----------------------//
//  Medium Coefficients  //
//----------------------//

/**
 * @brief Scalar diffusion coefficient in the heat equation
 */
class MediumScalarCoeff
        : public mfem::Coefficient
{
public:
    /**
     * @brief Constructor
     * @param testCase test case for the heat equation
     */
    MediumScalarCoeff(std::shared_ptr<TestCases> testCase)
        : m_testCase (testCase) {}

    //! Evaluates the scalar coefficient
    virtual double Eval(mfem::ElementTransformation &,
                        const mfem::IntegrationPoint &);

protected:
    std::shared_ptr<TestCases> m_testCase;
};


/**
 * @brief Matrix diffusion coefficient in the heat equation
 */
class MediumTensorCoeff
        : public mfem::MatrixCoefficient
{
public:
    /**
     * @brief Constructor
     * @param testCase test case for the heat equation
     */
    MediumTensorCoeff(std::shared_ptr<TestCases> testCase)
        : MatrixCoefficient (testCase->getDim()),
          m_testCase (testCase) {}

    //! Evaluates the matrix coefficient
    virtual void Eval(mfem::DenseMatrix &, mfem::ElementTransformation &,
                      const mfem::IntegrationPoint &);

protected:
    std::shared_ptr<TestCases> m_testCase;
};


//-------------------------------//
//  Exact Solution Coefficients  //
//-------------------------------//
/**
 * @brief Exact temperature (scalar) coefficient
 */
class ExactTemperatureCoeff
        : public mfem::Coefficient
{
public:
    /**
     * @brief Constructor
     * @param testCase test case for the heat equation
     */
    ExactTemperatureCoeff(std::shared_ptr<TestCases> testCase)
        : m_testCase (testCase) {}
    
    //! Evaluates the scalar coefficient, time t is deduced by GetTime()
    virtual double Eval(mfem::ElementTransformation &,
                        const mfem::IntegrationPoint &);
    
    //! Evaluates the scalar coefficient, time t is passed as an argument
    double Eval(mfem::ElementTransformation &,
                const mfem::IntegrationPoint &,
                double);
    
protected:
    std::shared_ptr<TestCases> m_testCase;
};


/**
 * @brief Exact heat flux (vector) coefficient
 */
class ExactHeatFluxCoeff
        : public mfem::VectorCoefficient
{
public:
    /**
     * @brief Constructor
     * @param testCase test case for the heat equation
     */
    ExactHeatFluxCoeff (std::shared_ptr<TestCases> testCase)
        : mfem::VectorCoefficient (testCase->getDim()),
          m_testCase (testCase) {}

    //! Evaluates the vector coefficient, time t is deduced by GetTime()
    virtual void Eval(mfem::Vector&, mfem::ElementTransformation&,
                      const mfem::IntegrationPoint&);
    
    //! Evaluates the vector coefficient, time t is passed as an argument
    void Eval(mfem::Vector&, mfem::ElementTransformation&,
              const mfem::IntegrationPoint&, double);
    
private:
    std::shared_ptr<TestCases> m_testCase;
};

/**
 * @brief Exact temperature spatial gradient (vector) coefficient
 */
class ExactTemperatureSpatialGradCoeff
        : public mfem::VectorCoefficient
{
public:
    /**
     * @brief Constructor
     * @param testCase test case for the heat equation
     */
    ExactTemperatureSpatialGradCoeff
    (std::shared_ptr<TestCases> testCase)
        : mfem::VectorCoefficient (testCase->getDim()),
          m_testCase (testCase) {}

    //! Evaluates the vector coefficient, time t is deduced by GetTime()
    virtual void Eval(mfem::Vector&, mfem::ElementTransformation&,
                      const mfem::IntegrationPoint&);

    //! Evaluates the vector coefficient, time t is passed as an argument
    void Eval(mfem::Vector&, mfem::ElementTransformation&,
              const mfem::IntegrationPoint&, double);

private:
    std::shared_ptr<TestCases> m_testCase;
};

/**
 * @brief Exact temperature gradient with respect to time (scalar) coefficient
 */
class ExactTemperatureTemporalGradCoeff
        : public mfem::Coefficient
{
public:
    /**
     * @brief Constructor
     * @param testCase test case for the heat equation
     */
    ExactTemperatureTemporalGradCoeff
    (std::shared_ptr<TestCases> testCase)
        : m_testCase (testCase) {}

    //! Evaluates the scalar coefficient, time t is deduced by GetTime()
    virtual double Eval(mfem::ElementTransformation &,
                        const mfem::IntegrationPoint &);

    //! Evaluates the scalar coefficient, time t is passed as an argument
    double Eval(mfem::ElementTransformation &,
                const mfem::IntegrationPoint &,
                double);

protected:
    std::shared_ptr<TestCases> m_testCase;
};

/**
 * @brief Laplacian (scalar) coefficient
 */
class LaplacianCoeff
        : public mfem::Coefficient
{
public:
    /**
     * @brief Constructor
     * @param testCase test case for the heat equation
     */
    LaplacianCoeff(std::shared_ptr<TestCases> testCase)
        : m_testCase (testCase) {}

    //! Evaluates the scalar coefficient, time t is deduced by GetTime()
    virtual double Eval(mfem::ElementTransformation &,
                        const mfem::IntegrationPoint &);

    //! Evaluates the scalar coefficient, time t is passed as an argument
    double Eval(mfem::ElementTransformation &,
                const mfem::IntegrationPoint &,
                double);

protected:
    std::shared_ptr<TestCases> m_testCase;
};


//---------------------------------//
//  Initial Solution Coefficients  //
//---------------------------------//

/**
 * @brief Initial temperature (scalar) coefficient
 */
class InitialTemperatureCoeff
        : public mfem::Coefficient
{
public:
    /**
     * @brief Constructor
     * @param testCase test case for the heat equation
     */
    InitialTemperatureCoeff(std::shared_ptr<TestCases> testCase)
        : m_testCase (testCase) {}

    //! Evaluates the scalar coefficient
    virtual double Eval(mfem::ElementTransformation &,
                        const mfem::IntegrationPoint &);

protected:
    std::shared_ptr<TestCases> m_testCase;
};


//------------------------------------//
//  Boundary Temperature Coefficient  //
//------------------------------------//

/**
 * @brief Boundary temperature (scalar) coefficient
 */
class BdryTemperatureCoeff
        : public mfem::Coefficient
{
public:
    /**
     * @brief Constructor
     * @param testCase test case for the heat equation
     */
    BdryTemperatureCoeff (std::shared_ptr<TestCases> testCase)
        : m_testCase (testCase) {}

    //! Evaluates the scalar coefficient, time t is deduced by GetTime()
    virtual double Eval(mfem::ElementTransformation &,
                        const mfem::IntegrationPoint &);

    //! Evaluates the scalar coefficient, time t is passed as an argument
    double Eval(mfem::ElementTransformation &,
                const mfem::IntegrationPoint &,
                double);

private:
    std::shared_ptr<TestCases> m_testCase;
};


//-----------------------//
//  Source Coefficients  //
//-----------------------//

/**
 * @brief Source (scalar) coefficient in the heat equation
 */
class SourceCoeff
        : public mfem::Coefficient
{
public:
    /**
     * @brief Constructor
     * @param testCase test case for the heat equation
     */
    SourceCoeff (std::shared_ptr<TestCases> testCase)
        : m_testCase (testCase) {}

    //! Evaluates the scalar coefficient, time t is deduced by GetTime()
    virtual double Eval(mfem::ElementTransformation &,
                        const mfem::IntegrationPoint &);

    //! Evaluates the scalar coefficient, time t is passed as an argument
    double Eval(mfem::ElementTransformation &,
                const mfem::IntegrationPoint &,
                double);

private:
    std::shared_ptr<TestCases> m_testCase;
};

}

#endif // HEAT_COEFFICIENTS_HPP
