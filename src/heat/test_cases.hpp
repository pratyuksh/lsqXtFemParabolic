#ifndef HEAT_TEST_CASES_HPP
#define HEAT_TEST_CASES_HPP

#include "../core/config.hpp"
#include "mfem.hpp"

using namespace mfem;

// Test cases available
enum {Dummy,
      UnitSquareTest1,
      UnitSquareTest2,
      UnitSquareTest3,
      UnitSquareTest4,
      PeriodicUnitSquareTest1,
      LShapedTest1,
      LShapedTest2};

namespace heat {

/**
 * @brief Abstract Base class for test cases of the heat equation
 */
class TestCases
{
public:
    //! Default constructor
    virtual ~TestCases() = default;
    
    /**
     * @brief Defines the scalar material coefficient
     * @return material coefficient at a given physical point
     */
    virtual double medium(const Vector&) const = 0;

    /**
     * @brief Defines the matrix material coefficient
     * @return material coefficient at a given physical point
     */
    virtual DenseMatrix mediumTensor(const Vector&) const = 0;

    /**
     * @brief Temperature solution
     * @return temperature value at a given point in space-time
     */
    virtual double temperatureSol(const Vector&,
                                  const double) const = 0;

    /**
     * @brief Heat flux solution
     * @return heat flux value at a given point in space-time
     */
    virtual Vector heatFluxSol(const Vector&,
                               const double) const = 0;

    /**
     * @brief Temperature gradient with respect to time
     * @return value of temperature time-gradient at a given point in space-time
     */
    virtual double temperatureTimeGradientSol (const Vector&,
                                               const double) const = 0;

    virtual double laplacian
    (const Vector&, const double) const = 0;
    
    /**
     * @brief Initial temperature
     * @return value of initial temperature at a given point in space
     */
    virtual double initTemperature(const Vector&) const = 0;

    /**
     * @brief Initial heat flux
     * @return value of heat flux vector at a given point in space
     */
    virtual Vector initHeatFlux(const Vector&) const = 0;
    
    /**
     * @brief Temperature at the boundary of spatial domain
     * @return value of temperature at a given point on the spatial boundary
     */
    virtual double bdryTemperature(const Vector&,
                                   const double) const = 0;

    /**
     * @brief source term in the heat equation
     * @return value of source at a given point in space-time
     */
    virtual double source(const Vector&,
                          const double) const = 0;
    
    //! Sets the Dirichlet boundary
    virtual void setBdryDirichlet(Array<int>&) const = 0;

    virtual void setPerturbation(double) const = 0;

    virtual void setPerturbations(const Vector&) const = 0;

    //! Returns the number of spatial dimensions
    virtual int getDim() const = 0;
};


/**
 * @brief Base template for different test cases
 */
template<int ProblemType>
class TestCase;


/**
 * Template specialization for a dummy test case;
 * value 1 everywhere, except boundary;
 * Homogeneous Dirichlet boundary;
 * No source/forcing
 */
template<>
class TestCase <Dummy>
        : public TestCases
{
public:
    explicit TestCase (const nlohmann::json& config)
        : m_config(config), m_dim(2) {
    }

    double temperatureSol(const Vector&,
                           const double) const override;

    Vector heatFluxSol(const Vector&,
                        const double) const override;

    double temperatureTimeGradientSol
        (const Vector&, const double) const override;

    double source(const Vector&,
                  const double) const override;

    void setBdryDirichlet(Array<int>&) const override;

    double medium(const Vector&) const override;

    DenseMatrix mediumTensor(const Vector&) const override;

    double laplacian
    (const Vector& x, const double t) const override {
        return temperatureTimeGradientSol(x,t) - source(x, t);
    }

    double initTemperature(const Vector& x) const override {
        return temperatureSol(x,0);
    }

    Vector initHeatFlux(const Vector& x) const override {
        return heatFluxSol(x,0);
    }

    double bdryTemperature(const Vector& x,
                            const double t) const override {
        return temperatureSol(x, t);
    }

    void setPerturbation(double) const override {}

    void setPerturbations(const Vector&) const override {}

    inline int getDim() const override {
        return m_dim;
    }

private:
    const nlohmann::json& m_config;
    int m_dim;
};


/**
 * Template specialization for a unit-square domain, with Test1;
 * Non-zero initial conditions;
 * Homogeneous Dirichlet boundary;
 * No source/forcing
 */
template<>
class TestCase <UnitSquareTest1>
        : public TestCases
{
public:
    explicit TestCase (const nlohmann::json& config)
        : m_config(config), m_dim(2) {
    }

    double temperatureSol(const Vector&,
                           const double) const override;

    Vector heatFluxSol(const Vector&,
                        const double) const override;

    double temperatureTimeGradientSol
        (const Vector&, const double) const override;

    double source(const Vector&,
                  const double) const override;

    void setBdryDirichlet(Array<int>&) const override;

    double medium(const Vector&) const override;

    DenseMatrix mediumTensor(const Vector&) const override;

    double laplacian
    (const Vector& x, const double t) const override {
        return temperatureTimeGradientSol(x,t) - source(x, t);
    }

    double initTemperature(const Vector& x) const override {
        return temperatureSol(x,0);
    }

    Vector initHeatFlux(const Vector& x) const override {
        return heatFluxSol(x,0);
    }

    double bdryTemperature(const Vector& x,
                            const double t) const override {
        return temperatureSol(x, t);
    }

    void setPerturbation(double w) const override {
        m_rvar = w;
    }

    void setPerturbations(const Vector&) const override {}

    inline int getDim() const override {
        return m_dim;
    }

private:
    const nlohmann::json& m_config;
    int m_dim;

    mutable double m_rvar = 0;
};

/**
 * Template specialization for a unit-square domain, with Test2;
 * Non-zero initial conditions;
 * Homogeneous Dirichlet boundary;
 * Non-zero forcing/source
 */
template<>
class TestCase <UnitSquareTest2>
        : public TestCases
{
public:
    explicit TestCase (const nlohmann::json& config)
        : m_config(config), m_dim(2) {
    }

    double temperatureSol(const Vector&,
                           const double) const override;

    Vector heatFluxSol(const Vector&,
                        const double) const override;

    double temperatureTimeGradientSol
        (const Vector&, const double) const override;

    double source(const Vector&,
                  const double) const override;

    void setBdryDirichlet(Array<int>&) const override;

    double medium(const Vector&) const override
    {
        return 1;
    }

    DenseMatrix mediumTensor(const Vector&) const override
    {
        DenseMatrix med(m_dim);
        med(0,0) = med(1,1) = 1;
        med(0,1) = med(1,0) = 0;
        return med;
    }

    double laplacian
    (const Vector& x, const double t) const override {
        return temperatureTimeGradientSol(x,t)
                - source(x, t);
    }

    double initTemperature(const Vector& x) const override {
        return temperatureSol(x,0);
    }

    Vector initHeatFlux(const Vector& x) const override {
        return heatFluxSol(x,0);
    }

    double bdryTemperature(const Vector& x,
                            const double t) const override {
        return temperatureSol(x, t);
    }

    void setPerturbation(double) const override {}

    void setPerturbations(const Vector&) const override {}

    inline int getDim() const override {
        return m_dim;
    }

private:
    const nlohmann::json& m_config;
    int m_dim;
};

/**
 * Template specialization for a unit-square domain, with Test3;
 * Zero initial conditions;
 * Homogeneous Dirichlet boundary;
 * Non-zero forcing/source
 */
template<>
class TestCase <UnitSquareTest3>
        : public TestCases
{
public:
    explicit TestCase (const nlohmann::json& config)
        : m_config(config), m_dim(2) {
    }

    double medium(const Vector&) const override;

    DenseMatrix mediumTensor(const Vector&) const override;

    double temperatureSol(const Vector&,
                           const double) const override;

    Vector heatFluxSol(const Vector&,
                        const double) const override;

    double temperatureTimeGradientSol
        (const Vector&, const double) const override;

    double source(const Vector&,
                  const double) const override;

    void setBdryDirichlet(Array<int>&) const override;

    double laplacian
    (const Vector& x, const double t) const override {
        return temperatureTimeGradientSol(x,t)
                - source(x, t);
    }

    double initTemperature(const Vector& x) const override {
        return temperatureSol(x,0);
    }

    Vector initHeatFlux(const Vector& x) const override {
        return heatFluxSol(x,0);
    }

    double bdryTemperature(const Vector& x,
                            const double t) const override {
        return temperatureSol(x, t);
    }

    void setPerturbation(double w) const override {
        m_rvar = w;
    }

    void setPerturbations(const Vector&) const override {}

    inline int getDim() const override {
        return m_dim;
    }

private:
    double perturb(const Vector&) const;

private:
    const nlohmann::json& m_config;
    int m_dim;

    mutable double m_rvar = 0;
};

/**
 * Template specialization for a unit-square domain, with Test4
 * Zero initial conditions;
 * Homogeneous Dirichlet boundary;
 * Non-zero forcing/source;
 * Material coefficient matrix
 */
template<>
class TestCase <UnitSquareTest4>
        : public TestCases
{
public:
    explicit TestCase (const nlohmann::json& config)
        : m_config(config), m_dim(2) {
    }

    DenseMatrix mediumTensor(const Vector&) const override;

    double temperatureSol(const Vector&,
                           const double) const override;

    Vector heatFluxSol(const Vector&,
                        const double) const override;

    double temperatureTimeGradientSol
        (const Vector&, const double) const override;

    double source(const Vector&,
                  const double) const override;

    void setBdryDirichlet(Array<int>&) const override;

    double laplacian
    (const Vector& x, const double t) const override {
        return temperatureTimeGradientSol(x,t)
                - source(x, t);
    }

    double initTemperature(const Vector& x) const override {
        return temperatureSol(x,0);
    }

    Vector initHeatFlux(const Vector& x) const override {
        return heatFluxSol(x,0);
    }

    double bdryTemperature(const Vector& x,
                            const double t) const override {
        return temperatureSol(x, t);
    }

    void setPerturbation(double) const override {}

    void setPerturbations(const Vector&) const override {}

    inline int getDim() const override {
        return m_dim;
    }

    // redundant for this test case
    double medium(const Vector&) const override
    {
        return 1;
    }

private:
    const nlohmann::json& m_config;
    int m_dim;
};

/**
 * Template specialization for a periodic unit-square domain, with Test1;
 * Zero initial conditions;
 * periodic boundary;
 * Non-zero source
 */
template<>
class TestCase <PeriodicUnitSquareTest1>
        : public TestCases
{
public:
    explicit TestCase (const nlohmann::json& config)
        : m_config(config), m_dim(2) {
    }

    double medium(const Vector& x) const override;

    DenseMatrix mediumTensor(const Vector&) const override;

    double temperatureSol(const Vector&,
                           const double) const override;

    Vector heatFluxSol(const Vector&,
                        const double) const override;

    double temperatureTimeGradientSol
        (const Vector&, const double) const override;

    double source(const Vector&,
                  const double) const override;

    void setPerturbation(double w) const override {
        m_rvar = w;
    }

    void setPerturbations(const Vector& w) const override {
        m_rvars = w;
    }

    double laplacian
    (const Vector& x, const double t) const override {
        return temperatureTimeGradientSol(x,t)
                - source(x, t);
    }

    double initTemperature(const Vector& x) const override {
        return temperatureSol(x,0);
    }

    Vector initHeatFlux(const Vector& x) const override {
        return heatFluxSol(x,0);
    }

    double bdryTemperature(const Vector& x,
                            const double t) const override {
        return temperatureSol(x, t);
    }

    void setBdryDirichlet(Array<int>&) const override {}

    inline int getDim() const override {
        return m_dim;
    }

private:
    double perturb(const Vector&) const;

private:
    const nlohmann::json& m_config;
    int m_dim;

    mutable double m_rvar = 0;
    mutable Vector m_rvars;
};

/**
 * Template specialization for an L-shaped domain, with Test1;
 * Non-zero initial conditions;
 * Homogeneous Dirichlet boundary;
 * Non-zero source
 */
template<>
class TestCase <LShapedTest1>
        : public TestCases
{
public:
    explicit TestCase (const nlohmann::json& config)
        : m_config(config), m_dim(2) {
        m_singularCorner.SetSize(m_dim);
        m_singularCorner = 0.0;
    }

    double temperatureSol(const Vector&,
                           const double) const override;

    Vector heatFluxSol(const Vector&,
                        const double) const override;

    double temperatureTimeGradientSol
        (const Vector&, const double) const override;

    double source(const Vector&,
                  const double) const override;

    void setBdryDirichlet(Array<int>&) const override;

    double medium(const Vector&) const override
    {
        return 1;
    }

    DenseMatrix mediumTensor(const Vector&) const override
    {
        DenseMatrix med(m_dim);
        med(0,0) = med(1,1) = 1;
        med(0,1) = med(1,0) = 0;
        return med;
    }

    double laplacian
    (const Vector& x, const double t) const override {
        return temperatureTimeGradientSol(x,t)
                - source(x, t);
    }

    double initTemperature(const Vector& x) const override {
        return temperatureSol(x,0);
    }

    Vector initHeatFlux(const Vector& x) const override {
        return heatFluxSol(x,0);
    }

    double bdryTemperature(const Vector& x,
                            const double t) const override {
        return temperatureSol(x, t);
    }

    void setPerturbation(double) const override {}

    void setPerturbations(const Vector&) const override {}

    inline int getDim() const override {
        return m_dim;
    }

private:
    double radius(const Vector& x) const;
    double polarAngle(const Vector& x) const;

    const nlohmann::json& m_config;
    int m_dim;

    double m_gamma = 2./3.;
    Vector m_singularCorner;
};

/**
 * Template specialization for an L-shaped domain, with Test2;
 * Zero initial conditions;
 * Homogeneous Dirichlet boundary;
 * Non-zero source
 */
template<>
class TestCase <LShapedTest2>
        : public TestCases
{
public:
    explicit TestCase (const nlohmann::json& config)
        : m_config(config), m_dim(2) {
        m_singularCorner.SetSize(m_dim);
        m_singularCorner = 0.0;
    }

    double temperatureSol(const Vector&,
                           const double) const override;

    Vector heatFluxSol(const Vector&,
                        const double) const override;

    double temperatureTimeGradientSol
        (const Vector&, const double) const override;

    double source(const Vector&,
                  const double) const override;

    void setBdryDirichlet(Array<int>&) const override;

    double medium(const Vector&) const override
    {
        return 1;
    }

    DenseMatrix mediumTensor(const Vector&) const override
    {
        DenseMatrix med(m_dim);
        med(0,0) = med(1,1) = 1;
        med(0,1) = med(1,0) = 0;
        return med;
    }

    double laplacian
    (const Vector& x, const double t) const override {
        return temperatureTimeGradientSol(x,t)
                - source(x, t);
    }

    double initTemperature(const Vector& x) const override {
        return temperatureSol(x,0);
    }

    Vector initHeatFlux(const Vector& x) const override {
        return heatFluxSol(x,0);
    }

    double bdryTemperature(const Vector& x,
                            const double t) const override {
        return temperatureSol(x, t);
    }

    void setPerturbation(double) const override {}

    void setPerturbations(const Vector&) const override {}

    inline int getDim() const override {
        return m_dim;
    }

private:
    double radius(const Vector& x) const;
    double polarAngle(const Vector& x) const;

    const nlohmann::json& m_config;
    int m_dim;

    double m_gamma = 2./3.;
    Vector m_singularCorner;
};

}

#endif // HEAT_TEST_CASES_HPP
