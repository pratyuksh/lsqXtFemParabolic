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
class MediumScalarCoeff : public Coefficient
{
public:
    /**
     * @brief Constructor
     * @param testCase test case for the heat equation
     */
    MediumScalarCoeff(std::shared_ptr<TestCases> testCase)
        : m_testCase (testCase) {}

    //! Evaluates the scalar coefficient
    virtual double Eval(ElementTransformation &,
                        const IntegrationPoint &);

protected:
    std::shared_ptr<TestCases> m_testCase;
};


/**
 * @brief Matrix diffusion coefficient in the heat equation
 */
class MediumTensorCoeff : public MatrixCoefficient
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
    virtual void Eval(DenseMatrix &, ElementTransformation &,
                      const IntegrationPoint &);

protected:
    std::shared_ptr<TestCases> m_testCase;
};


//-------------------------------//
//  Exact Solution Coefficients  //
//-------------------------------//
/**
 * @brief Exact temperature (scalar) coefficient
 */
class ExactTemperatureCoeff : public Coefficient
{
public:
    /**
     * @brief Constructor
     * @param testCase test case for the heat equation
     */
    ExactTemperatureCoeff(std::shared_ptr
                              <TestCases> testCase)
        : m_testCase (testCase) {}
    
    //! Evaluates the scalar coefficient, time t is deduced by GetTime()
    virtual double Eval(ElementTransformation &,
                        const IntegrationPoint &);
    
    //! Evaluates the scalar coefficient, time t is passed as an argument
    double Eval(ElementTransformation &,
                const IntegrationPoint &,
                double);
    
protected:
    std::shared_ptr<TestCases> m_testCase;
};


/**
 * @brief Exact heat flux (vector) coefficient
 */
class ExactFluxCoeff : public VectorCoefficient
{
public:
    /**
     * @brief Constructor
     * @param testCase test case for the heat equation
     */
    ExactFluxCoeff (std::shared_ptr
                        <TestCases> testCase)
        : VectorCoefficient (testCase->getDim()),
          m_testCase (testCase) {}

    //! Evaluates the vector coefficient, time t is deduced by GetTime()
    virtual void Eval(Vector&, ElementTransformation&,
                      const IntegrationPoint&);
    
    //! Evaluates the vector coefficient, time t is passed as an argument
    void Eval(Vector&, ElementTransformation&,
              const IntegrationPoint&, double);
    
private:
    std::shared_ptr<TestCases> m_testCase;
};

/**
 * @brief Exact temperature gradient with respect to time (scalar) coefficient
 */
class ExactTemperatureTimeGradCoeff : public Coefficient
{
public:
    /**
     * @brief Constructor
     * @param testCase test case for the heat equation
     */
    ExactTemperatureTimeGradCoeff
    (std::shared_ptr<TestCases> testCase)
        : m_testCase (testCase) {}

    //! Evaluates the scalar coefficient, time t is deduced by GetTime()
    virtual double Eval(ElementTransformation &,
                        const IntegrationPoint &);

    //! Evaluates the scalar coefficient, time t is passed as an argument
    double Eval(ElementTransformation &,
                const IntegrationPoint &,
                double);

protected:
    std::shared_ptr<TestCases> m_testCase;
};

/**
 * @brief Laplacian (scalar) coefficient
 */
class LaplacianCoeff : public Coefficient
{
public:
    /**
     * @brief Constructor
     * @param testCase test case for the heat equation
     */
    LaplacianCoeff
    (std::shared_ptr<TestCases> testCase)
        : m_testCase (testCase) {}

    //! Evaluates the scalar coefficient, time t is deduced by GetTime()
    virtual double Eval(ElementTransformation &,
                        const IntegrationPoint &);

    //! Evaluates the scalar coefficient, time t is passed as an argument
    double Eval(ElementTransformation &,
                const IntegrationPoint &,
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
class InitialTemperatureCoeff : public Coefficient
{
public:
    /**
     * @brief Constructor
     * @param testCase test case for the heat equation
     */
    InitialTemperatureCoeff(std::shared_ptr
                              <TestCases> testCase)
        : m_testCase (testCase) {}

    //! Evaluates the scalar coefficient
    virtual double Eval(ElementTransformation &,
                        const IntegrationPoint &);

protected:
    std::shared_ptr<TestCases> m_testCase;
};


//------------------------------------//
//  Boundary Temperature Coefficient  //
//------------------------------------//

/**
 * @brief Boundary temperature (scalar) coefficient
 */
class BdryTemperatureCoeff : public Coefficient
{
public:
    /**
     * @brief Constructor
     * @param testCase test case for the heat equation
     */
    BdryTemperatureCoeff (std::shared_ptr
                              <TestCases> testCase)
        : m_testCase (testCase) {}

    //! Evaluates the scalar coefficient, time t is deduced by GetTime()
    virtual double Eval(ElementTransformation &,
                        const IntegrationPoint &);

    //! Evaluates the scalar coefficient, time t is passed as an argument
    double Eval(ElementTransformation &,
                const IntegrationPoint &,
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
class SourceCoeff : public Coefficient
{
public:
    /**
     * @brief Constructor
     * @param testCase test case for the heat equation
     */
    SourceCoeff (std::shared_ptr
                     <TestCases> testCase)
        : m_testCase (testCase) {}

    //! Evaluates the scalar coefficient, time t is deduced by GetTime()
    virtual double Eval(ElementTransformation &,
                        const IntegrationPoint &);

    //! Evaluates the scalar coefficient, time t is passed as an argument
    double Eval(ElementTransformation &,
                const IntegrationPoint &,
                double);

private:
    std::shared_ptr<TestCases> m_testCase;
};

}

#endif // HEAT_COEFFICIENTS_HPP
