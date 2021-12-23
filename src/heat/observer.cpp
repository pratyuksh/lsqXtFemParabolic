#include "observer.hpp"

#include <fstream>
#include <iostream>
#include <filesystem>

#include "../mymfem/utilities.hpp"
#include "utilities.hpp"

using namespace mfem;
namespace fs = std::filesystem;


// Constructor
heat::Observer
:: Observer (const nlohmann::json& config, int spatialLevel)
    : BaseObserver (config, spatialLevel)
{
    m_endTime = 1;
    if (config.contains("end_time")) {
        m_endTime = config["end_time"];
    }

    m_boolEvalError = false;
    if (config.contains("eval_error")) {
        m_boolEvalError = config["eval_error"];
    }

    m_errorType = "natural";
    if (config.contains("error_type")) {
        m_errorType = config["error_type"];
    }
    
    std::string problem_type = config["problem_type"];
    m_outputDir = "../output/heat_"+problem_type+"/";
    if (m_boolDumpOut) {
        fs::create_directories(m_outputDir);
    }
    m_solNameSuffix = "_lx"+std::to_string(spatialLevel);
}

// Sets the test case, discretisation, etc.
void heat::Observer
:: set (std::shared_ptr<heat::TestCases>& testCase,
        std::shared_ptr<heat::LsqXtFem>& disc)
{
    m_testCase = testCase;
    m_disc = disc;

    m_xDim = m_disc->getSpatialMesh()->Dimension();
    m_temporalFeSpace = m_disc->getTemporalFeSpace();
    m_spatialFeSpaceForTemperature = m_disc->getSpatialFeSpaces()[0];
    m_spatialFeSpaceForHeatFlux = m_disc->getSpatialFeSpaces()[1];

    if (m_boolEvalError)
    {
        m_materialCoeff = std::make_unique<heat::MediumTensorCoeff>
                (m_testCase);

        m_exactTemperatureCoeff = std::make_unique
                <heat::ExactTemperatureCoeff>(m_testCase);
        m_exactHeatFluxCoeff = std::make_unique<heat::ExactHeatFluxCoeff>
                (m_testCase);

        m_exactSpatialGradientOfTemperatureCoeff = std::make_unique
                <heat::ExactTemperatureSpatialGradCoeff>
                (m_testCase);
        m_exactTemporalGradientOfTemperatureCoeff = std::make_unique
                <heat::ExactTemperatureTemporalGradCoeff>
                (m_testCase);

        m_sourceCoeff = std::make_unique<heat::SourceCoeff>
                (m_testCase);

        m_exactTemperature = std::make_unique<GridFunction>
                (m_spatialFeSpaceForTemperature);
        m_exactHeatFlux = std::make_unique<GridFunction>
                (m_spatialFeSpaceForHeatFlux);
    }
}

void heat::Observer
:: visualizeSolutionAtEndTime
(const heat::SolutionHandler& solutionHandler)
{
    auto spatialFeSpaces = m_disc->getSpatialFeSpaces();

    auto temperatureDataAtEndTime
            = solutionHandler.getTemperatureDataAtEndTime();

    auto temperatureAtEndTime
            = std::make_unique<GridFunction>
            (spatialFeSpaces[0],
            temperatureDataAtEndTime.GetData());
    visualize(*temperatureAtEndTime);

    auto heatFluxDataAtEndTime
            = solutionHandler.getHeatFluxDataAtEndTime();

    auto heatFluxAtEndTime
            = std::make_unique<GridFunction>
            (spatialFeSpaces[1],
            heatFluxDataAtEndTime.GetData());
    visualize(*heatFluxAtEndTime);
}

// Writes the solution to output file
void heat::Observer
:: dumpSol (std::shared_ptr<GridFunction>& temperature,
            std::shared_ptr<GridFunction>& heatFlux) const
{
    if (m_boolDumpOut)
    {
        std::string solName1
                = m_outputDir+"temperature"
                +m_solNameSuffix;
        std::string solName2
                = m_outputDir+"heatFlux"+m_solNameSuffix;
        std::cout << solName1 << std::endl;
        std::cout << solName2 << std::endl;

        std::ofstream solOfs1(solName1.c_str());
        solOfs1.precision(m_precision);
        temperature->Save(solOfs1);

        std::ofstream solOfs2(solName2.c_str());
        solOfs2.precision(m_precision);
        heatFlux->Save(solOfs2);
    }
}

Vector heat::Observer
:: evalError(const heat::SolutionHandler& solutionHandler)
{
    Vector solutionError;

    if (m_errorType == "natural") {
        solutionError = evalErrorInNaturalNorm(solutionHandler);
    }
    else if (m_errorType == "lsq") {
        solutionError = evalErrorInLeastSquaresNorm(solutionHandler);
    }

    return solutionError;
}


Vector heat::Observer
:: evalErrorInNaturalNorm
(const heat::SolutionHandler& solutionHandler) const
{
    Vector solutionError(3);

    auto temperatureData = solutionHandler.getTemperatureData();
    auto heatFluxData = solutionHandler.getHeatFluxData();

    Vector temperatureRelErrorNormL2H1(2);
    Vector heatFluxRelErrorNormL2L2(2);
    Vector solDivergenceRelErrorNormL2L2(2);

    int numTemporalMeshEls = m_temporalFeSpace->GetNE();
    int spatialFeSpaceSizeForTemperature
            = m_spatialFeSpaceForTemperature->GetTrueVSize();
    int spatialFeSpaceSizeForHeatFlux
            = m_spatialFeSpaceForHeatFlux->GetTrueVSize();

    Vector temperatureSol(spatialFeSpaceSizeForTemperature);
    std::shared_ptr<GridFunction> temperatureSolGF
            = std::make_shared<GridFunction>
            (m_spatialFeSpaceForTemperature, temperatureSol);

    Vector temporalGradientOfTemperatureSol(spatialFeSpaceSizeForTemperature);
    std::shared_ptr<GridFunction> temporalGradientOfTemperatureSolGF
            = std::make_shared<GridFunction>
            (m_spatialFeSpaceForTemperature, temporalGradientOfTemperatureSol);

    Vector heatFluxSol(spatialFeSpaceSizeForHeatFlux);
    std::shared_ptr<GridFunction> heatFluxSolGF
            = std::make_shared<GridFunction>
            (m_spatialFeSpaceForHeatFlux, heatFluxSol);

    ElementTransformation *temporalTrans = nullptr;
    const FiniteElement *temporalFe = nullptr;
    Vector temporalShape;
    DenseMatrix temporalDShape;
    Vector temporalGrad;
    Array<int> temporalVdofs;

    // auxiliary variables
    double localTemperatureErrorNormH1, localTemperatureNormH1;
    double localHeatFluxErrorNormL2, localHeatFluxNormL2;
    double localSolDivergenceErrorNormL2, localSolDivergenceNormL2;

    temperatureRelErrorNormL2H1 = 0.;
    heatFluxRelErrorNormL2L2 = 0.;
    solDivergenceRelErrorNormL2L2 = 0.;
    for (int n=0; n<numTemporalMeshEls; n++)
    {
        m_temporalFeSpace->GetElementVDofs(n, temporalVdofs);
        temporalTrans = m_temporalFeSpace->GetElementTransformation(n);

        // compute elvec
        temporalFe = m_temporalFeSpace->GetFE(n);
        int tNdofs = temporalFe->GetDof();
        temporalShape.SetSize(tNdofs);
        temporalDShape.SetSize(tNdofs, 1);
        temporalDShape.GetColumnReference(0, temporalGrad);

        int order = 2*temporalFe->GetOrder()+1;
        const IntegrationRule *ir
                = &IntRules.Get(temporalFe->GetGeomType(), order);

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir->IntPoint(i);
            temporalTrans->SetIntPoint(&ip);
            temporalFe->CalcShape(ip, temporalShape);
            temporalFe->CalcPhysDShape(*temporalTrans, temporalDShape);

            // build solution at time t
            Vector t;
            temporalTrans->Transform(ip, t);
            buildSolutionAtSpecifiedTime(temperatureData, temporalShape,
                                         temporalVdofs, temperatureSol);
            buildSolutionAtSpecifiedTime(temperatureData, temporalGrad,
                                         temporalVdofs, temporalGradientOfTemperatureSol);
            buildSolutionAtSpecifiedTime(heatFluxData, temporalShape,
                                         temporalVdofs, heatFluxSol);

            std::tie (localTemperatureErrorNormH1, localTemperatureNormH1,
                      localHeatFluxErrorNormL2, localHeatFluxNormL2,
                      localSolDivergenceErrorNormL2, localSolDivergenceNormL2)
                    =  evalSpatialErrorOfSolutionInNaturalNorm
                    (temperatureSolGF, temporalGradientOfTemperatureSolGF,
                     heatFluxSolGF, t(0));

            double w = ip.weight*temporalTrans->Weight();

            temperatureRelErrorNormL2H1(0)
                    += w*localTemperatureErrorNormH1*localTemperatureErrorNormH1;
            temperatureRelErrorNormL2H1(1)
                    += w*localTemperatureNormH1*localTemperatureNormH1;

            heatFluxRelErrorNormL2L2(0)
                    += w*localHeatFluxErrorNormL2*localHeatFluxErrorNormL2;
            heatFluxRelErrorNormL2L2(1)
                    += w*localHeatFluxNormL2*localHeatFluxNormL2;

            solDivergenceRelErrorNormL2L2(0)
                    += w*localSolDivergenceErrorNormL2*localSolDivergenceErrorNormL2;
            solDivergenceRelErrorNormL2L2(1)
                    += w*localSolDivergenceNormL2*localSolDivergenceNormL2;
        }
    }
    solutionError(0) = std::sqrt(temperatureRelErrorNormL2H1(0))
            / std::sqrt(temperatureRelErrorNormL2H1(1));
    solutionError(1) = std::sqrt(heatFluxRelErrorNormL2L2(0))
            / std::sqrt(heatFluxRelErrorNormL2L2(1));
    solutionError(2) = std::sqrt(solDivergenceRelErrorNormL2L2(0))
            / std::sqrt(solDivergenceRelErrorNormL2L2(1));

    return solutionError;
}

// Evaluates error in the natural norm at a given time t
std::tuple <double, double, double, double, double, double>
heat::Observer
:: evalSpatialErrorOfSolutionInNaturalNorm
(std::shared_ptr<GridFunction>& temperature,
 std::shared_ptr<GridFunction>& temporalGradientOfTemperature,
 std::shared_ptr<GridFunction>& heatFlux,
 double t) const
{
    double temperatureErrorNormH1=0, temperatureNormH1=0;
    double heatFluxErrorNormL2=0, heatFluxNormL2=0;
    double solDivergenceNormErrorL2=0, solDivergenceNormL2=0;

    auto spatialMesh = temperature->FESpace()->GetMesh();
    m_exactTemperatureCoeff->SetTime(t);
    m_exactHeatFluxCoeff->SetTime(t);
    m_exactSpatialGradientOfTemperatureCoeff->SetTime(t);
    m_exactTemporalGradientOfTemperatureCoeff->SetTime(t);
    m_sourceCoeff->SetTime(t);

    // loop over all elements of the space mesh,
    Vector temperatureSol;
    DenseMatrix spatialGradientOfTemperatureSol;
    DenseMatrix heatFluxSol;
    const FiniteElement *fe = nullptr;
    ElementTransformation *trans= nullptr;
    for (int i=0; i<spatialMesh->GetNE(); i++)
    {
        fe = m_spatialFeSpaceForTemperature->GetFE(i);
        int order = 2*fe->GetOrder()+1;
        const IntegrationRule *ir
                = &IntRules.Get(fe->GetGeomType(), order);

        trans = m_spatialFeSpaceForTemperature->GetElementTransformation(i);
        temperature->GetValues(i, *ir, temperatureSol);
        temperature->GetGradients(*trans, *ir,
                                  spatialGradientOfTemperatureSol);
        heatFlux->GetVectorValues(*trans, *ir, heatFluxSol);

        Vector localSpatialGradientOfTemperatureErrorNormL2;
        Vector localHeatFluxErrorNormL2;
        int numPoints = ir->GetNPoints();
        for (int j=0; j<numPoints; j++)
        {
            const IntegrationPoint &ip = ir->IntPoint(j);
            trans->SetIntPoint(&ip);
            double w = ip.weight*trans->Weight();

            // H1-error temperature
            double exactTemperatureSol
                    = m_exactTemperatureCoeff->Eval(*trans, ip);
            m_exactSpatialGradientOfTemperatureCoeff
                    ->Eval(localSpatialGradientOfTemperatureErrorNormL2,
                           *trans, ip);
            temperatureNormH1
                    += w*(exactTemperatureSol*exactTemperatureSol
                          + (localSpatialGradientOfTemperatureErrorNormL2
                             *localSpatialGradientOfTemperatureErrorNormL2));
            double localTemperatureErrorNormL2
                    = exactTemperatureSol - temperatureSol(j);
            Vector tmp1(spatialGradientOfTemperatureSol.GetColumn(j),
                        m_xDim);
            localSpatialGradientOfTemperatureErrorNormL2.Add(-1, tmp1);
            temperatureErrorNormH1
                    += w*(localTemperatureErrorNormL2*localTemperatureErrorNormL2 +
                          localSpatialGradientOfTemperatureErrorNormL2
                          *localSpatialGradientOfTemperatureErrorNormL2);

            // L2-error heat-flux
            m_exactHeatFluxCoeff->Eval(localHeatFluxErrorNormL2, *trans, ip);
            heatFluxNormL2 += w*(localHeatFluxErrorNormL2
                                 *localHeatFluxErrorNormL2);
            Vector tmp2(heatFluxSol.GetColumn(j), m_xDim);
            localHeatFluxErrorNormL2.Add(-1, tmp2);
            heatFluxErrorNormL2 += w*(localHeatFluxErrorNormL2
                                      *localHeatFluxErrorNormL2);

            // solution divergence error
            double localSolDivergenceNormL2
                    = m_sourceCoeff->Eval(*trans, ip);
            solDivergenceNormL2 += w*(localSolDivergenceNormL2
                                      *localSolDivergenceNormL2);
            localSolDivergenceNormL2 -=
                    (temporalGradientOfTemperature->GetValue(i, ip)
                     - heatFlux->GetDivergence(*trans));
            solDivergenceNormErrorL2 += w*(localSolDivergenceNormL2
                                           *localSolDivergenceNormL2);
        }
    }
    temperatureErrorNormH1 = std::sqrt(temperatureErrorNormH1);
    temperatureNormH1 = std::sqrt(temperatureNormH1);

    heatFluxErrorNormL2 = std::sqrt(heatFluxErrorNormL2);
    heatFluxNormL2 = std::sqrt(heatFluxNormL2);

    solDivergenceNormErrorL2 = std::sqrt(solDivergenceNormErrorL2);
    solDivergenceNormL2 = std::sqrt(solDivergenceNormL2);

    return {std::move(temperatureErrorNormH1), std::move(temperatureNormH1),
                std::move(heatFluxErrorNormL2), std::move(heatFluxNormL2),
                std::move(solDivergenceNormErrorL2), std::move(solDivergenceNormL2)};
}

Vector heat::Observer
:: evalErrorInLeastSquaresNorm
(const heat::SolutionHandler& solutionHandler) const
{
    Vector solutionError(3);

    auto temperatureData = solutionHandler.getTemperatureData();
    auto heatFluxData = solutionHandler.getHeatFluxData();

    int numTemporalMeshEls = m_temporalFeSpace->GetNE();
    int spatialFeSpaceSizeForTemperature
            = m_spatialFeSpaceForTemperature->GetTrueVSize();
    int spatialFeSpaceSizeForHeatFlux
            = m_spatialFeSpaceForHeatFlux->GetTrueVSize();

    Vector temperatureSol(spatialFeSpaceSizeForTemperature);
    std::shared_ptr<GridFunction> temperatureSolGF
            = std::make_shared<GridFunction>
            (m_spatialFeSpaceForTemperature, temperatureSol);

    Vector temporalGradientOfTemperatureSol(spatialFeSpaceSizeForTemperature);
    std::shared_ptr<GridFunction> temporalGradientOfTemperatureSolGF
            = std::make_shared<GridFunction>
            (m_spatialFeSpaceForTemperature, temporalGradientOfTemperatureSol);

    Vector heatFluxSol(spatialFeSpaceSizeForHeatFlux);
    std::shared_ptr<GridFunction> heatFluxSolGF
            = std::make_shared<GridFunction>
            (m_spatialFeSpaceForHeatFlux, heatFluxSol);

    ElementTransformation *temporalTrans = nullptr;
    const FiniteElement *temporalFe = nullptr;
    Vector temporalShape;
    DenseMatrix temporalDShape;
    Vector temporalGrad;
    Array<int> temporalVdofs;

    // error in PDE and flux
    double pdeErrorNormL2L2 = 0;
    double heatFluxErrorNormL2L2 = 0;

    double localPdeErrorNormL2, localHeatFluxErrorNormL2;
    for (int n=0; n<numTemporalMeshEls; n++)
    {
        m_temporalFeSpace->GetElementVDofs(n, temporalVdofs);
        temporalTrans = m_temporalFeSpace->GetElementTransformation(n);

        // compute elvec
        temporalFe = m_temporalFeSpace->GetFE(n);
        int tNdofs = temporalFe->GetDof();

        temporalShape.SetSize(tNdofs);
        temporalDShape.SetSize(tNdofs, 1);
        temporalDShape.GetColumnReference(0, temporalGrad);

        int order = 2*temporalFe->GetOrder()+1;
        const IntegrationRule *ir
                = &IntRules.Get(temporalFe->GetGeomType(), order);

        for (int i = 0; i < ir->GetNPoints(); i++)
        {
            const IntegrationPoint &ip = ir->IntPoint(i);
            temporalTrans->SetIntPoint(&ip);
            temporalFe->CalcShape(ip, temporalShape);
            temporalFe->CalcPhysDShape(*temporalTrans, temporalDShape);

            // build solution at time t
            Vector t;
            temporalTrans->Transform(ip, t);
            buildSolutionAtSpecifiedTime(temperatureData, temporalShape,
                                         temporalVdofs, temperatureSol);
            buildSolutionAtSpecifiedTime(temperatureData, temporalGrad,
                                         temporalVdofs, temporalGradientOfTemperatureSol);
            buildSolutionAtSpecifiedTime(heatFluxData, temporalShape,
                                         temporalVdofs, heatFluxSol);

            std::tie (localPdeErrorNormL2, localHeatFluxErrorNormL2)
                    = evalSpatialErrorOfSolutionInLeastSquaresNorm
                    (temperatureSolGF, temporalGradientOfTemperatureSolGF,
                     heatFluxSolGF, t(0));

            double w = ip.weight*temporalTrans->Weight();
            pdeErrorNormL2L2 += w*localPdeErrorNormL2*localPdeErrorNormL2;
            heatFluxErrorNormL2L2 += w*localHeatFluxErrorNormL2
                    *localHeatFluxErrorNormL2;
        }
    }
    solutionError(0) = std::sqrt(pdeErrorNormL2L2);
    solutionError(1) = std::sqrt(heatFluxErrorNormL2L2);

    // error in initial conditions
    double initialTemperatureErrorNormL2 = 0;
    {
        int n=0;
        m_temporalFeSpace->GetElementVDofs(n, temporalVdofs);
        temporalTrans = m_temporalFeSpace->GetElementTransformation(n);

        temporalFe = m_temporalFeSpace->GetFE(n);
        int temporalNumDofs = temporalFe->GetDof();
        temporalShape.SetSize(temporalNumDofs);

        IntegrationPoint ip;
        ip.Set1w(0, 1.0);
        temporalFe->CalcShape(ip, temporalShape);

        buildSolutionAtSpecifiedTime(temperatureData, temporalShape,
                                     temporalVdofs, temperatureSol);

        m_exactTemperatureCoeff->SetTime(0);
        initialTemperatureErrorNormL2 = temperatureSolGF
                ->ComputeL2Error(*m_exactTemperatureCoeff);
    }
    solutionError(2) = initialTemperatureErrorNormL2;

    return solutionError;
}

std::tuple <double, double> heat::Observer
:: evalSpatialErrorOfSolutionInLeastSquaresNorm
(std::shared_ptr<GridFunction>& temperature,
 std::shared_ptr<GridFunction>& temporalGradientOfTemperature,
 std::shared_ptr<GridFunction>& heatFlux,
 double t) const
{
    double pdeErrorNormL2=0, heatFluxErrorNormL2=0;

    auto spatialMesh = temperature->FESpace()->GetMesh();
    m_sourceCoeff->SetTime(t);

    // loop over all elements of the space mesh,
    DenseMatrix gradTemperatureSol, heatFluxSol;
    DenseMatrix materialTensor;
    const FiniteElement *fe = nullptr;
    ElementTransformation *trans= nullptr;
    for (int i=0; i<spatialMesh->GetNE(); i++)
    {
        fe = m_spatialFeSpaceForTemperature->GetFE(i);
        int order = 2*fe->GetOrder()+1;
        const IntegrationRule *ir
                = &IntRules.Get(fe->GetGeomType(), order);

        trans = m_spatialFeSpaceForTemperature->GetElementTransformation(i);
        temperature->GetGradients(*trans, *ir, gradTemperatureSol);
        heatFlux->GetVectorValues(*trans, *ir, heatFluxSol);

        Vector localHeatFluxErrorNormL2;
        int numPoints = ir->GetNPoints();
        for (int j=0; j<numPoints; j++)
        {
            const IntegrationPoint &ip = ir->IntPoint(j);
            trans->SetIntPoint(&ip);
            double w = ip.weight*trans->Weight();

            // pde
            double localPdeErrorNormL2
                    = temporalGradientOfTemperature->GetValue(i, ip);
            localPdeErrorNormL2 -= heatFlux->GetDivergence(*trans);
            localPdeErrorNormL2 -= m_sourceCoeff->Eval(*trans, ip);
            pdeErrorNormL2 += w*(localPdeErrorNormL2*localPdeErrorNormL2);

            // flux
            m_materialCoeff->Eval(materialTensor, *trans, ip);
            heatFluxSol.GetColumnReference(j, localHeatFluxErrorNormL2);
            Vector tmp(gradTemperatureSol.GetColumn(j), m_xDim);
            materialTensor.AddMult_a(-1, tmp, localHeatFluxErrorNormL2);
            heatFluxErrorNormL2 += w*(localHeatFluxErrorNormL2
                                      *localHeatFluxErrorNormL2);
        }
    }
    pdeErrorNormL2 = std::sqrt(pdeErrorNormL2);
    heatFluxErrorNormL2 = std::sqrt(heatFluxErrorNormL2);

    return {pdeErrorNormL2, heatFluxErrorNormL2};
}

// End of file
