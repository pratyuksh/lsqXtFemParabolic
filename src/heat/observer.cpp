#include "observer.hpp"

#include <fstream>
#include <iostream>
#include <filesystem>

#include <Eigen/Core>

#include "../mymfem/utilities.hpp"
#include "utilities.hpp"
#include "functionals.hpp"

#include <omp.h>

namespace fs = std::filesystem;


// Constructor
heat::Observer
:: Observer (const nlohmann::json& config, int lx)
    : BaseObserver (config, lx)
{
    m_boolError = false;
    if (config.contains("eval_error")) {
        m_boolError = config["eval_error"];
    }
    
    std::string problem_type = config["problem_type"];
    m_outputDir = "../output/heat_"+problem_type+"/";
    if (m_boolDumpOut) {
        fs::create_directories(m_outputDir);
    }
    m_solNameSuffix = "_lx"+std::to_string(lx);
}

// Sets the test case, discretisation, etc.
void heat::Observer
:: set (std::shared_ptr<heat::TestCases>& testCase,
        std::shared_ptr<heat::LsqXtFEM>& discr)
{
    m_testCase = testCase;
    m_discr = discr;

    m_xndim = m_discr->getXMesh()->Dimension();
    m_tVspace = m_discr->getTFespace();
    m_xV1space = m_discr->getXFespaces()[0];
    m_xV2space = m_discr->getXFespaces()[1];

    if (m_boolError)
    {
        m_medCoeff = std::make_unique<heat::MediumTensorCoeff>
                (m_testCase);
        m_uECoeff = std::make_unique
                <heat::ExactTemperatureCoeff>(m_testCase);
        m_qECoeff = std::make_unique<heat::ExactFluxCoeff>
                (m_testCase);
        m_dudtECoeff = std::make_unique
                <heat::ExactTemperatureTimeGradCoeff>
                (m_testCase);
        m_sourceCoeff = std::make_unique<heat::SourceCoeff>
                (m_testCase);
        //m_laplacianCoeff = std::make_unique
        //        <heat::LaplacianCoeff> (m_testCase);

        m_uE = std::make_unique<GridFunction>(m_xV1space);
        m_qE = std::make_unique<GridFunction>(m_xV2space);
    }

    // observational and quantity of interest
    m_obsFn = std::make_unique<QoIFunctional>(m_xV1space);
    m_qoiFn = std::make_unique<QoIFunctional>
            (m_xV1space, 1.5);

    m_qoiXtFn = std::make_unique<QoIXtFunctional>
            (m_tVspace, m_xV1space, 1.5);
}

// Dumps solution
void heat::Observer
:: dumpSol (std::shared_ptr<GridFunction>& u,
             std::shared_ptr<GridFunction>& q) const
{
    if (m_boolDumpOut)
    {
        std::string sol_name1
                = m_outputDir+"temperature"
                +m_solNameSuffix;
        std::string sol_name2
                = m_outputDir+"heatFlux"+m_solNameSuffix;
        std::cout << sol_name1 << std::endl;
        std::cout << sol_name2 << std::endl;

        std::ofstream sol_ofs1(sol_name1.c_str());
        sol_ofs1.precision(m_precision);
        u->Save(sol_ofs1);

        std::ofstream sol_ofs2(sol_name2.c_str());
        sol_ofs2.precision(m_precision);
        q->Save(sol_ofs2);
    }
}

// Computes L2 and H1 error for temperature
std::tuple <double, double, double, double> heat::Observer
:: evalXH1Error (std::shared_ptr<GridFunction>& u,
                  double t) const
{
    double erruL2=0, uEL2=0;
    double erruH1=0, uEH1=0;

    ConstantCoefficient zero(0);
    ConstantCoefficient one(1.0);
    VectorFunctionCoefficient zeroVec(m_xndim, zeroFn);

    // L2 error in temperature
    m_uECoeff->SetTime(t);
    m_uE->ProjectCoefficient(*m_uECoeff);
    uEL2 = m_uE->ComputeL2Error(zero);
    erruL2 = u->ComputeL2Error(*m_uECoeff);

    // H1 error in temperature
    m_qECoeff->SetTime(t);
    uEH1 = m_uE->ComputeH1Error(&zero, &zeroVec, &one,
                                1.0, 1);
    erruH1 = u->ComputeH1Error(m_uECoeff.get(),
                               m_qECoeff.get(),
                               &one, 1.0, 1);

    return {std::move(erruL2), std::move(uEL2),
                std::move(erruH1), std::move(uEH1)};
}

// Computes L2L2 and L2H1 error for temperature
std::tuple <double, double, double, double> heat::Observer
:: evalXtL2H1Error (Vector& u) const
{
    double erruL2L2=0, uEL2L2=0;
    double erruL2H1=0, uEL2H1=0;

    int Nt = m_tVspace->GetNE();
    int xdimV = m_xV1space->GetTrueVSize();

    Vector uSol(xdimV);
    std::shared_ptr<GridFunction> uGSol
            = std::make_shared<GridFunction>
            (m_xV1space, uSol);

    ElementTransformation *tTrans = nullptr;
    const FiniteElement *tFe = nullptr;
    Vector tShape;

    Eigen::VectorXd bufErruL2L2(Nt); bufErruL2L2.setZero();
    Eigen::VectorXd bufErruL2H1(Nt); bufErruL2H1.setZero();
    Eigen::VectorXd bufuEL2L2(Nt); bufuEL2L2.setZero();
    Eigen::VectorXd bufuEL2H1(Nt); bufuEL2H1.setZero();

    Array<int> tVdofs;
    double erruL2, uEL2, erruH1, uEH1;

    for (int n=0; n<Nt; n++)
    {
        m_tVspace->GetElementVDofs(n, tVdofs);
        tTrans = m_tVspace->GetElementTransformation(n);

        // compute elvec
        tFe = m_tVspace->GetFE(n);
        int tNdofs = tFe->GetDof();
        tShape.SetSize(tNdofs);

        int order = 2*tFe->GetOrder()+1;
        const IntegrationRule *ir
                = &IntRules.Get(tFe->GetGeomType(), order);

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir->IntPoint(i);
            tTrans->SetIntPoint(&ip);
            tFe->CalcShape(ip, tShape);

            // build solution at time t
            Vector t;
            tTrans->Transform(ip, t);
            buildXSolFG(u, tShape, tVdofs, uSol);

            std::tie (erruL2, uEL2, erruH1, uEH1)
                    = evalXH1Error(uGSol, t(0));

            double w = ip.weight*tTrans->Weight();
            bufErruL2L2(n) += w*erruL2*erruL2;
            bufErruL2H1(n) += w*erruH1*erruH1;
            bufuEL2L2(n) += w*uEL2*uEL2;
            bufuEL2H1(n) += w*uEH1*uEH1;
        }
    }
    erruL2L2 = std::sqrt(bufErruL2L2.sum());
    erruL2H1 = std::sqrt(bufErruL2H1.sum());
    uEL2L2 = std::sqrt(bufuEL2L2.sum());
    uEL2H1 = std::sqrt(bufuEL2H1.sum());

    bufErruL2L2.setZero();
    bufErruL2H1.setZero();
    bufuEL2L2.setZero();
    bufuEL2H1.setZero();

    return {std::move(erruL2L2), std::move(uEL2L2),
                std::move(erruL2H1), std::move(uEL2H1)};
}

//! Computes least-squares error components
//! given a time t
std::tuple <double, double> heat::Observer
:: evalXLsqError (std::shared_ptr<GridFunction>& u,
                   std::shared_ptr<GridFunction>& dudt,
                   std::shared_ptr<GridFunction>& q,
                   double t) const
{
    double errPde=0, errFlux=0;

    auto xMesh = u->FESpace()->GetMesh();
    m_sourceCoeff->SetTime(t);

    // loop over all elements of the space mesh,
    DenseMatrix graduxSol, fluxSol;
    DenseMatrix medTensor;
    const FiniteElement *fe = nullptr;
    ElementTransformation *trans= nullptr;
    for (int i=0; i<xMesh->GetNE(); i++)
    {
        fe = m_xV1space->GetFE(i);
        int order = 2*fe->GetOrder()+1;
        const IntegrationRule *ir
                = &IntRules.Get(fe->GetGeomType(), order);

        trans = m_xV1space->GetElementTransformation(i);
        u->GetGradients(*trans, *ir, graduxSol);
        q->GetVectorValues(*trans, *ir, fluxSol);

        Vector errFluxLocal;
        int numPoints = ir->GetNPoints();
        for (int j=0; j<numPoints; j++)
        {
            const IntegrationPoint &ip = ir->IntPoint(j);
            trans->SetIntPoint(&ip);
            double w = ip.weight*trans->Weight();

            // pde
            double errPdeLocal = dudt->GetValue(i, ip);
            errPdeLocal -= q->GetDivergence(*trans);
            errPdeLocal -= m_sourceCoeff->Eval(*trans, ip);
            errPde += w*(errPdeLocal*errPdeLocal);

            // flux
            m_medCoeff->Eval(medTensor, *trans, ip);
            fluxSol.GetColumnReference(j, errFluxLocal);
            Vector tmp(graduxSol.GetColumn(j), m_xndim);
            medTensor.AddMult_a(-1, tmp, errFluxLocal);
            errFlux += w*(errFluxLocal*errFluxLocal);
        }
    }
    errPde = std::sqrt(errPde);
    errFlux = std::sqrt(errFlux);

    return {std::move(errPde), std::move(errFlux)};
}

//! Computes least-squares error components
std::tuple <double, double, double> heat::Observer
:: evalXtLsqError (Vector& u, Vector& q) const
{
    double errPde=0, errFlux=0, errIc=0;

    int Nt = m_tVspace->GetNE();
    int xdimV1 = m_xV1space->GetTrueVSize();
    int xdimV2 = m_xV2space->GetTrueVSize();

    Vector uSol(xdimV1);
    std::shared_ptr<GridFunction> uGSol
            = std::make_shared<GridFunction>
            (m_xV1space, uSol);

    Vector dudtSol(xdimV1);
    std::shared_ptr<GridFunction> dudtGSol
            = std::make_shared<GridFunction>
            (m_xV1space, dudtSol);

    Vector qSol(xdimV2);
    std::shared_ptr<GridFunction> qGSol
            = std::make_shared<GridFunction>
            (m_xV2space, qSol);

    ElementTransformation *tTrans = nullptr;
    const FiniteElement *tFe = nullptr;
    Vector tShape;
    Array<int> tVdofs;

    // error in initial conditions
    {
        int n=0;
        m_tVspace->GetElementVDofs(n, tVdofs);
        tTrans = m_tVspace->GetElementTransformation(n);

        tFe = m_tVspace->GetFE(n);
        int tNdofs = tFe->GetDof();
        tShape.SetSize(tNdofs);

        IntegrationPoint ip;
        ip.Set1w(0, 1.0);
        tFe->CalcShape(ip, tShape);

        buildXSolFG(u, tShape, tVdofs, uSol);
        //(*this)(uGSol);

        m_uECoeff->SetTime(0);
        errIc = uGSol->ComputeL2Error(*m_uECoeff);
    }

    // error in PDE and flux
    Eigen::VectorXd bufErrPde(Nt); bufErrPde.setZero();
    Eigen::VectorXd bufErrFlux(Nt); bufErrFlux.setZero();
    double errPdeLocal=0, errFluxLocal=0;
    DenseMatrix tDShape;
    Vector tGrad;
    for (int n=0; n<Nt; n++)
    {
        m_tVspace->GetElementVDofs(n, tVdofs);
        tTrans = m_tVspace->GetElementTransformation(n);

        // compute elvec
        tFe = m_tVspace->GetFE(n);
        int tNdofs = tFe->GetDof();
        tShape.SetSize(tNdofs);
        tDShape.SetSize(tNdofs, 1);
        tDShape.GetColumnReference(0, tGrad);

        int order = 2*tFe->GetOrder()+1;
        const IntegrationRule *ir
                = &IntRules.Get(tFe->GetGeomType(), order);

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir->IntPoint(i);
            tTrans->SetIntPoint(&ip);
            tFe->CalcShape(ip, tShape);
            tFe->CalcPhysDShape(*tTrans, tDShape);

            // build solution at time t
            Vector t;
            tTrans->Transform(ip, t);
            buildXSolFG(u, tShape, tVdofs, uSol);
            buildXSolFG(u, tGrad, tVdofs, dudtSol);
            buildXSolFG(q, tShape, tVdofs, qSol);

            std::tie (errPdeLocal, errFluxLocal)
                    = evalXLsqError(uGSol, dudtGSol,
                                     qGSol, t(0));

            double w = ip.weight*tTrans->Weight();
            bufErrPde(n) += w*errPdeLocal*errPdeLocal;
            bufErrFlux(n) += w*errFluxLocal*errFluxLocal;
        }
    }
    errPde = std::sqrt(bufErrPde.sum());
    errFlux = std::sqrt(bufErrFlux.sum());

    return {std::move(errPde), std::move(errFlux),
                std::move(errIc)};
}

//! Computes the error in the natural norm
//! at a given time t
std::tuple <double, double, double, double, double, double>
heat::Observer
:: evalXError (std::shared_ptr<GridFunction>& u,
                std::shared_ptr<GridFunction>& dudt,
                std::shared_ptr<GridFunction>& q,
                double t) const
{
    double erruH1=0, uEH1=0;
    double errqL2=0, qEL2=0;
    double errUDiv=0, UEDiv=0;

    auto xMesh = u->FESpace()->GetMesh();
    m_uECoeff->SetTime(t);
    m_qECoeff->SetTime(t);
    m_dudtECoeff->SetTime(t);
    m_sourceCoeff->SetTime(t);

    // H1-error in u
    ConstantCoefficient zero(0);
    ConstantCoefficient one(1.0);
    VectorFunctionCoefficient zeroVec(m_xndim, zeroFn);

    m_uE->ProjectCoefficient(*m_uECoeff);
    uEH1 = m_uE->ComputeH1Error(&zero, &zeroVec, &one,
                                1.0, 1);
    erruH1 = u->ComputeH1Error(m_uECoeff.get(),
                               m_qECoeff.get(),
                               &one, 1.0, 1);

    // loop over all elements of the space mesh,
    DenseMatrix graduxSol, fluxSol;
    DenseMatrix medTensor;
    const FiniteElement *fe = nullptr;
    ElementTransformation *trans= nullptr;
    for (int i=0; i<xMesh->GetNE(); i++)
    {
        fe = m_xV1space->GetFE(i);
        int order = 2*fe->GetOrder()+1;
        const IntegrationRule *ir
                = &IntRules.Get(fe->GetGeomType(), order);

        trans = m_xV1space->GetElementTransformation(i);
        u->GetGradients(*trans, *ir, graduxSol);
        q->GetVectorValues(*trans, *ir, fluxSol);

        Vector errqL2Local;
        int numPoints = ir->GetNPoints();
        for (int j=0; j<numPoints; j++)
        {
            const IntegrationPoint &ip = ir->IntPoint(j);
            trans->SetIntPoint(&ip);
            double w = ip.weight*trans->Weight();

            // L2-error q
            m_qECoeff->Eval(errqL2Local, *trans, ip);
            Vector tmp(fluxSol.GetColumn(j), m_xndim);
            qEL2 += w*(errqL2Local*errqL2Local);
            errqL2Local.Add(-1, tmp);
            errqL2 += w*(errqL2Local*errqL2Local);

            // divergence error
            double errUDivLocal
                    = m_sourceCoeff->Eval(*trans, ip);
            UEDiv += w*(errUDivLocal*errUDivLocal);
            errUDivLocal -= (dudt->GetValue(i, ip)
                             -q->GetDivergence(*trans));
            errUDiv += w*(errUDivLocal*errUDivLocal);
        }
    }
    errqL2 = std::sqrt(errqL2);   qEL2 = std::sqrt(qEL2);
    errUDiv = std::sqrt(errUDiv); UEDiv = std::sqrt(UEDiv);

    return {std::move(erruH1), std::move(uEH1),
                std::move(errqL2), std::move(qEL2),
                std::move(errUDiv), std::move(UEDiv)};
}

//! Computes the error in the natural norm
std::tuple <double, double, double> heat::Observer
:: evalXtError (Vector& u, Vector& q) const
{
    double erruL2H1=0, uEL2H1=0;
    double errqL2L2=0, qEL2L2=0;
    double errUDiv=0, UEDiv=0;

    int Nt = m_tVspace->GetNE();
    int xdimV1 = m_xV1space->GetTrueVSize();
    int xdimV2 = m_xV2space->GetTrueVSize();

    Vector uSol(xdimV1);
    std::shared_ptr<GridFunction> uGSol
            = std::make_shared<GridFunction>
            (m_xV1space, uSol);

    Vector dudtSol(xdimV1);
    std::shared_ptr<GridFunction> dudtGSol
            = std::make_shared<GridFunction>
            (m_xV1space, dudtSol);

    Vector qSol(xdimV2);
    std::shared_ptr<GridFunction> qGSol
            = std::make_shared<GridFunction>
            (m_xV2space, qSol);

    ElementTransformation *tTrans = nullptr;
    const FiniteElement *tFe = nullptr;
    Vector tShape;
    Array<int> tVdofs;

    // auxiliary variables
    Eigen::VectorXd bufErruL2H1(Nt); bufErruL2H1.setZero();
    Eigen::VectorXd bufuEL2H1(Nt); bufuEL2H1.setZero();

    Eigen::VectorXd bufErrqL2L2(Nt); bufErrqL2L2.setZero();
    Eigen::VectorXd bufqEL2L2(Nt); bufqEL2L2.setZero();

    Eigen::VectorXd bufErrUDiv(Nt); bufErrUDiv.setZero();
    Eigen::VectorXd bufUEDiv(Nt); bufUEDiv.setZero();

    double erruH1=0, uEH1=0;
    double errqL2=0, qEL2=0;
    double errUDivLocal=0, UEDivLocal=0;
    DenseMatrix tDShape;
    Vector tGrad;
    for (int n=0; n<Nt; n++)
    {
        m_tVspace->GetElementVDofs(n, tVdofs);
        tTrans = m_tVspace->GetElementTransformation(n);

        // compute elvec
        tFe = m_tVspace->GetFE(n);
        int tNdofs = tFe->GetDof();
        tShape.SetSize(tNdofs);
        tDShape.SetSize(tNdofs, 1);
        tDShape.GetColumnReference(0, tGrad);

        int order = 2*tFe->GetOrder()+1;
        const IntegrationRule *ir
                = &IntRules.Get(tFe->GetGeomType(), order);

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir->IntPoint(i);
            tTrans->SetIntPoint(&ip);
            tFe->CalcShape(ip, tShape);
            tFe->CalcPhysDShape(*tTrans, tDShape);

            // build solution at time t
            Vector t;
            tTrans->Transform(ip, t);
            buildXSolFG(u, tShape, tVdofs, uSol);
            buildXSolFG(u, tGrad, tVdofs, dudtSol);
            buildXSolFG(q, tShape, tVdofs, qSol);

            std::tie (erruH1, uEH1, errqL2, qEL2,
                      errUDivLocal, UEDivLocal)
                    =  evalXError(uGSol, dudtGSol,
                                   qGSol, t(0));

            double w = ip.weight*tTrans->Weight();

            bufErruL2H1(n) += w*erruH1*erruH1;
            bufuEL2H1(n) += w*uEH1*uEH1;

            bufErrqL2L2(n) += w*errqL2*errqL2;
            bufqEL2L2(n) += w*qEL2*qEL2;

            bufErrUDiv(n) += w*errUDivLocal*errUDivLocal;
            bufUEDiv(n) += w*UEDivLocal*UEDivLocal;
        }
    }
    erruL2H1 = std::sqrt(bufErruL2H1.sum());
    uEL2H1 = std::sqrt(bufuEL2H1.sum());
    erruL2H1 /= uEL2H1;

    errqL2L2 = std::sqrt(bufErrqL2L2.sum());
    qEL2L2 = std::sqrt(bufqEL2L2.sum());
    errqL2L2 /= qEL2L2;

    errUDiv = std::sqrt(bufErrUDiv.sum());
    UEDiv = std::sqrt(bufUEDiv.sum());
    errUDiv /= UEDiv;

    return {std::move(erruL2H1), std::move(errqL2L2),
                std::move(errUDiv)};
}

// Computes the flux observation functional at end time T
// given the temperature solution U in space-time
double heat::Observer
:: evalObs(Vector& U) const
{
    // dofs for the time element (T-h, T)
    Array<int> tVdofs;
    const FiniteElement *tFe
            = m_tVspace->GetFE(m_tVspace->GetNE()-1);
    m_tVspace->GetElementVDofs
            (m_tVspace->GetNE()-1, tVdofs);

    Vector tShape(tVdofs.Size());
    IntegrationPoint ip;
    ip.Set1w(1.0, 1.0);
    tFe->CalcShape(ip, tShape);

    // flux solution at T
    int xdimV = m_xV1space->GetTrueVSize();
    Vector uSol(xdimV);
    buildXSolFG(U, tShape, tVdofs, uSol);

    return std::move((*m_obsFn)(uSol));
}

// Computes the flux observation functional
// given the temperature solution at any time t
double heat::Observer
:: evalObs(std::shared_ptr <GridFunction>& u) const
{
    return std::move((*m_obsFn)(*u));
}

// Computes the quantity of interest functional at end time T
// given the temperature solution U in space-time
double heat::Observer
:: evalQoi(Vector& U) const
{
    // dofs for the time element (T-h, T)
    Array<int> tVdofs;
    const FiniteElement *tFe
            = m_tVspace->GetFE(m_tVspace->GetNE()-1);
    m_tVspace->GetElementVDofs
            (m_tVspace->GetNE()-1, tVdofs);

    Vector tShape(tVdofs.Size());
    IntegrationPoint ip;
    ip.Set1w(1.0, 1.0);
    tFe->CalcShape(ip, tShape);

    // temperature solution at T
    int xdimV = m_xV1space->GetTrueVSize();
    Vector uSol(xdimV);
    buildXSolFG(U, tShape, tVdofs, uSol);

    return std::move((*m_qoiFn)(uSol));
}

// Computes the quantity of interest functional
// given the temperature solution at any time t
double heat::Observer
:: evalQoi(std::shared_ptr <GridFunction>& u) const
{
    return std::move((*m_qoiFn)(*u));
}

double heat::Observer
:: evalXtQoi(Vector& U) const
{
    return std::move((*m_qoiXtFn)(U));
}


// End of file
