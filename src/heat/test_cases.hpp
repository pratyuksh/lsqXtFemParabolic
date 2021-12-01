#ifndef HEAT_TEST_CASES_HPP
#define HEAT_TEST_CASES_HPP

#include "../core/config.hpp"
#include "mfem.hpp"

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
    virtual double medium(const mfem::Vector&) const = 0;

    /**
     * @brief Defines the matrix material coefficient
     * @return material coefficient at a given physical point
     */
    virtual mfem::DenseMatrix mediumTensor(const mfem::Vector&) const = 0;

    /**
     * @brief Temperature solution
     * @return temperature value at a given point in space-time
     */
    virtual double temperatureSol(const mfem::Vector&,
                                  const double) const = 0;

    /**
     * @brief Heat flux solution
     * @return heat flux value at a given point in space-time
     */
    virtual mfem::Vector heatFluxSol(const mfem::Vector&,
                                     const double) const = 0;

    /**
     * @brief Temperature gradient with respect to time
     * @return value of temperature time-gradient at a given point in space-time
     */
    virtual double temperatureTimeGradientSol (const mfem::Vector&,
                                               const double) const = 0;

    virtual double laplacian
    (const mfem::Vector&, const double) const = 0;
    
    /**
     * @brief Initial temperature
     * @return value of initial temperature at a given point in space
     */
    virtual double initTemperature(const mfem::Vector&) const = 0;

    /**
     * @brief Initial heat flux
     * @return value of heat flux vector at a given point in space
     */
    virtual mfem::Vector initHeatFlux(const mfem::Vector&) const = 0;
    
    /**
     * @brief Temperature at the boundary of spatial domain
     * @return value of temperature at a given point on the spatial boundary
     */
    virtual double bdryTemperature(const mfem::Vector&,
                                   const double) const = 0;

    /**
     * @brief source term in the heat equation
     * @return value of source at a given point in space-time
     */
    virtual double source(const mfem::Vector&,
                          const double) const = 0;
    
    //! Sets the Dirichlet boundary
    virtual void setBdryDirichlet(mfem::Array<int>&) const = 0;

    virtual void setPerturbation(double) const = 0;

    virtual void setPerturbations(const mfem::Vector&) const = 0;

    //! Returns the number of spatial dimensions
    int getDim() const {
        return m_dim;
    }

protected:
    int m_dim;
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
        : m_config(config) {
        m_dim = 2;
    }

    double temperatureSol(const mfem::Vector&,
                          const double) const override;

    mfem::Vector heatFluxSol(const mfem::Vector&,
                             const double) const override;

    double temperatureTimeGradientSol
        (const mfem::Vector&, const double) const override;

    double source(const mfem::Vector&,
                  const double) const override;

    void setBdryDirichlet(mfem::Array<int>&) const override;

    double medium(const mfem::Vector&) const override;

    mfem::DenseMatrix mediumTensor(const mfem::Vector&) const override;

    double laplacian
    (const mfem::Vector& x, const double t) const override {
        return temperatureTimeGradientSol(x,t) - source(x, t);
    }

    double initTemperature(const mfem::Vector& x) const override {
        return temperatureSol(x,0);
    }

    mfem::Vector initHeatFlux(const mfem::Vector& x) const override {
        return heatFluxSol(x,0);
    }

    double bdryTemperature(const mfem::Vector& x,
                            const double t) const override {
        return temperatureSol(x, t);
    }

    void setPerturbation(double) const override {}

    void setPerturbations(const mfem::Vector&) const override {}

private:
    const nlohmann::json& m_config;
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
        : m_config(config) {
        m_dim = 2;
    }

    double temperatureSol(const mfem::Vector&,
                          const double) const override;

    mfem::Vector heatFluxSol(const mfem::Vector&,
                             const double) const override;

    double temperatureTimeGradientSol
        (const mfem::Vector&, const double) const override;

    double source(const mfem::Vector&,
                  const double) const override;

    void setBdryDirichlet(mfem::Array<int>&) const override;

    double medium(const mfem::Vector&) const override;

    mfem::DenseMatrix mediumTensor(const mfem::Vector&) const override;

    double laplacian
    (const mfem::Vector& x, const double t) const override {
        return temperatureTimeGradientSol(x,t) - source(x, t);
    }

    double initTemperature(const mfem::Vector& x) const override {
        return temperatureSol(x,0);
    }

    mfem::Vector initHeatFlux(const mfem::Vector& x) const override {
        return heatFluxSol(x,0);
    }

    double bdryTemperature(const mfem::Vector& x,
                           const double t) const override {
        return temperatureSol(x, t);
    }

    void setPerturbation(double w) const override {
        m_rvar = w;
    }

    void setPerturbations(const mfem::Vector&) const override {}

private:
    const nlohmann::json& m_config;
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
        : m_config(config) {
        m_dim = 2;
    }

    double temperatureSol(const mfem::Vector&,
                          const double) const override;

    mfem::Vector heatFluxSol(const mfem::Vector&,
                             const double) const override;

    double temperatureTimeGradientSol
        (const mfem::Vector&, const double) const override;

    double source(const mfem::Vector&,
                  const double) const override;

    void setBdryDirichlet(mfem::Array<int>&) const override;

    double medium(const mfem::Vector&) const override {
        return 1;
    }

    mfem::DenseMatrix mediumTensor(const mfem::Vector&) const override
    {
        mfem::DenseMatrix med(m_dim);
        med(0,0) = med(1,1) = 1;
        med(0,1) = med(1,0) = 0;
        return med;
    }

    double laplacian
    (const mfem::Vector& x, const double t) const override {
        return temperatureTimeGradientSol(x,t)
                - source(x, t);
    }

    double initTemperature(const mfem::Vector& x) const override {
        return temperatureSol(x,0);
    }

    mfem::Vector initHeatFlux(const mfem::Vector& x) const override {
        return heatFluxSol(x,0);
    }

    double bdryTemperature(const mfem::Vector& x,
                           const double t) const override {
        return temperatureSol(x, t);
    }

    void setPerturbation(double) const override {}

    void setPerturbations(const mfem::Vector&) const override {}

private:
    const nlohmann::json& m_config;
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
        : m_config(config) {
        m_dim = 2;
    }

    double medium(const mfem::Vector&) const override;

    mfem::DenseMatrix mediumTensor(const mfem::Vector&) const override;

    double temperatureSol(const mfem::Vector&,
                          const double) const override;

    mfem::Vector heatFluxSol(const mfem::Vector&,
                             const double) const override;

    double temperatureTimeGradientSol
        (const mfem::Vector&, const double) const override;

    double source(const mfem::Vector&,
                  const double) const override;

    void setBdryDirichlet(mfem::Array<int>&) const override;

    double laplacian
    (const mfem::Vector& x, const double t) const override {
        return temperatureTimeGradientSol(x,t)
                - source(x, t);
    }

    double initTemperature(const mfem::Vector& x) const override {
        return temperatureSol(x,0);
    }

    mfem::Vector initHeatFlux(const mfem::Vector& x) const override {
        return heatFluxSol(x,0);
    }

    double bdryTemperature(const mfem::Vector& x,
                           const double t) const override {
        return temperatureSol(x, t);
    }

    void setPerturbation(double w) const override {
        m_rvar = w;
    }

    void setPerturbations(const mfem::Vector&) const override {}

private:
    double perturb(const mfem::Vector&) const;

private:
    const nlohmann::json& m_config;
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
        : m_config(config) {
        m_dim = 2;
    }

    mfem::DenseMatrix mediumTensor(const mfem::Vector&) const override;

    double temperatureSol(const mfem::Vector&,
                          const double) const override;

    mfem::Vector heatFluxSol(const mfem::Vector&,
                             const double) const override;

    double temperatureTimeGradientSol
        (const mfem::Vector&, const double) const override;

    double source(const mfem::Vector&,
                  const double) const override;

    void setBdryDirichlet(mfem::Array<int>&) const override;

    double laplacian
    (const mfem::Vector& x, const double t) const override {
        return temperatureTimeGradientSol(x,t)
                - source(x, t);
    }

    double initTemperature(const mfem::Vector& x) const override {
        return temperatureSol(x,0);
    }

    mfem::Vector initHeatFlux(const mfem::Vector& x) const override {
        return heatFluxSol(x,0);
    }

    double bdryTemperature(const mfem::Vector& x,
                           const double t) const override {
        return temperatureSol(x, t);
    }

    void setPerturbation(double) const override {}

    void setPerturbations(const mfem::Vector&) const override {}

    // redundant for this test case
    double medium(const mfem::Vector&) const override {
        return 1;
    }

private:
    const nlohmann::json& m_config;
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
        : m_config(config) {
        m_dim = 2;
    }

    double medium(const mfem::Vector& x) const override;

    mfem::DenseMatrix mediumTensor(const mfem::Vector&) const override;

    double temperatureSol(const mfem::Vector&,
                          const double) const override;

    mfem::Vector heatFluxSol(const mfem::Vector&,
                             const double) const override;

    double temperatureTimeGradientSol
        (const mfem::Vector&, const double) const override;

    double source(const mfem::Vector&,
                  const double) const override;

    void setPerturbation(double w) const override {
        m_rvar = w;
    }

    void setPerturbations(const mfem::Vector& w) const override {
        m_rvars = w;
    }

    double laplacian
    (const mfem::Vector& x, const double t) const override {
        return temperatureTimeGradientSol(x,t)
                - source(x, t);
    }

    double initTemperature(const mfem::Vector& x) const override {
        return temperatureSol(x,0);
    }

    mfem::Vector initHeatFlux(const mfem::Vector& x) const override {
        return heatFluxSol(x,0);
    }

    double bdryTemperature(const mfem::Vector& x,
                           const double t) const override {
        return temperatureSol(x, t);
    }

    void setBdryDirichlet(mfem::Array<int>&) const override {}

private:
    double perturb(const mfem::Vector&) const;

private:
    const nlohmann::json& m_config;

    mutable double m_rvar = 0;
    mutable mfem::Vector m_rvars;
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
        : m_config(config) {
        m_dim = 2;
        m_singularCorner.SetSize(m_dim);
        m_singularCorner = 0.0;
    }

    double temperatureSol(const mfem::Vector&,
                          const double) const override;

    mfem::Vector heatFluxSol(const mfem::Vector&,
                             const double) const override;

    double temperatureTimeGradientSol
        (const mfem::Vector&, const double) const override;

    double source(const mfem::Vector&,
                  const double) const override;

    void setBdryDirichlet(mfem::Array<int>&) const override;

    double medium(const mfem::Vector&) const override {
        return 1;
    }

    mfem::DenseMatrix mediumTensor(const mfem::Vector&) const override
    {
        mfem::DenseMatrix med(m_dim);
        med(0,0) = med(1,1) = 1;
        med(0,1) = med(1,0) = 0;
        return med;
    }

    double laplacian
    (const mfem::Vector& x, const double t) const override {
        return temperatureTimeGradientSol(x,t)
                - source(x, t);
    }

    double initTemperature(const mfem::Vector& x) const override {
        return temperatureSol(x,0);
    }

    mfem::Vector initHeatFlux(const mfem::Vector& x) const override {
        return heatFluxSol(x,0);
    }

    double bdryTemperature(const mfem::Vector& x,
                           const double t) const override {
        return temperatureSol(x, t);
    }

    void setPerturbation(double) const override {}

    void setPerturbations(const mfem::Vector&) const override {}

private:
    double radius(const mfem::Vector& x) const;
    double polarAngle(const mfem::Vector& x) const;

private:
    const nlohmann::json& m_config;

    double m_gamma = 2./3.;
    mfem::Vector m_singularCorner;
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
        : m_config(config) {
        m_dim = 2;
        m_singularCorner.SetSize(m_dim);
        m_singularCorner = 0.0;
    }

    double temperatureSol(const mfem::Vector&,
                          const double) const override;

    mfem::Vector heatFluxSol(const mfem::Vector&,
                             const double) const override;

    double temperatureTimeGradientSol
        (const mfem::Vector&, const double) const override;

    double source(const mfem::Vector&,
                  const double) const override;

    void setBdryDirichlet(mfem::Array<int>&) const override;

    double medium(const mfem::Vector&) const override {
        return 1;
    }

    mfem::DenseMatrix mediumTensor(const mfem::Vector&) const override
    {
        mfem::DenseMatrix med(m_dim);
        med(0,0) = med(1,1) = 1;
        med(0,1) = med(1,0) = 0;
        return med;
    }

    double laplacian
    (const mfem::Vector& x, const double t) const override {
        return temperatureTimeGradientSol(x,t)
                - source(x, t);
    }

    double initTemperature(const mfem::Vector& x) const override {
        return temperatureSol(x,0);
    }

    mfem::Vector initHeatFlux(const mfem::Vector& x) const override {
        return heatFluxSol(x,0);
    }

    double bdryTemperature(const mfem::Vector& x,
                           const double t) const override {
        return temperatureSol(x, t);
    }

    void setPerturbation(double) const override {}

    void setPerturbations(const mfem::Vector&) const override {}

private:
    double radius(const mfem::Vector& x) const;
    double polarAngle(const mfem::Vector& x) const;

private:
    const nlohmann::json& m_config;

    double m_gamma = 2./3.;
    mfem::Vector m_singularCorner;
};

}

#endif // HEAT_TEST_CASES_HPP
