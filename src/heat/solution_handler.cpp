#include "solution_handler.hpp"

#include "utilities.hpp"

using namespace mfem;


heat::SolutionHandler
:: SolutionHandler
(mfem::FiniteElementSpace* temporalFeSpace,
 mfem::Array<mfem::FiniteElementSpace *> spatialFeSpaces)
    : m_temporalFeSpace (temporalFeSpace),
      m_spatialFeSpaces (spatialFeSpaces)
{
    m_temporalFeSpaceSize = temporalFeSpace->GetTrueVSize();
    m_spatialFeSpaceSizeForTemperature
            = spatialFeSpaces[0]->GetTrueVSize();
    m_spatialFeSpaceSizeForHeatFlux
            = spatialFeSpaces[1]->GetTrueVSize();

    int temperatureDataSize
            = m_temporalFeSpaceSize*m_spatialFeSpaceSizeForTemperature;
    int heatFluxDataSize
            = m_temporalFeSpaceSize*m_spatialFeSpaceSizeForHeatFlux;

    Array<int> blockOffsets(3);
    blockOffsets[0] = 0;
    blockOffsets[1] = temperatureDataSize;
    blockOffsets[2] = blockOffsets[1] + heatFluxDataSize;

    m_data = std::make_shared<BlockVector>(blockOffsets);
    (*m_data) = 0.;
}

Vector heat::SolutionHandler
:: getTemperatureDataAtInitialTime() const
{
    auto temperatureData = getTemperatureData();

    Array<int> temporalVdofs;
    const FiniteElement *temporalFe
            = m_temporalFeSpace->GetFE(0);
    m_temporalFeSpace->GetElementVDofs(0, temporalVdofs);

    Vector temporalShape(temporalVdofs.Size());
    IntegrationPoint ip;
    ip.Set1w(0.0, 1.0);
    temporalFe->CalcShape(ip, temporalShape);

    Vector temperatureAtInitialTime(m_spatialFeSpaceSizeForTemperature);
    buildSolutionAtSpecifiedTime(temperatureData,
                                 temporalShape,
                                 temporalVdofs,
                                 temperatureAtInitialTime);

    return temperatureAtInitialTime;
}

Vector heat::SolutionHandler
:: getTemperatureDataAtEndTime() const
{
    auto temporalMesh = m_temporalFeSpace->GetMesh();
    auto temperatureData = getTemperatureData();

    Array<int> temporalVdofs;
    const FiniteElement *temporalFe
            = m_temporalFeSpace->GetFE(temporalMesh->GetNE()-1);
    m_temporalFeSpace->GetElementVDofs(temporalMesh->GetNE()-1,
                                       temporalVdofs);

    Vector temporalShape(temporalVdofs.Size());
    IntegrationPoint ip;
    ip.Set1w(1.0, 1.0);
    temporalFe->CalcShape(ip, temporalShape);

    Vector temperatureAtEndTime(m_spatialFeSpaceSizeForTemperature);
    buildSolutionAtSpecifiedTime(temperatureData,
                                 temporalShape,
                                 temporalVdofs,
                                 temperatureAtEndTime);

    return temperatureAtEndTime;
}

Vector heat::SolutionHandler
:: getHeatFluxDataAtInitialTime() const
{
    auto heatFluxData = getHeatFluxData();

    Array<int> temporalVdofs;
    const FiniteElement *temporalFe
            = m_temporalFeSpace->GetFE(0);
    m_temporalFeSpace->GetElementVDofs(0, temporalVdofs);

    Vector temporalShape(temporalVdofs.Size());
    IntegrationPoint ip;
    ip.Set1w(0.0, 1.0);
    temporalFe->CalcShape(ip, temporalShape);

    Vector heatFluxAtInitialTime(m_spatialFeSpaceSizeForHeatFlux);
    buildSolutionAtSpecifiedTime(heatFluxData,
                                 temporalShape,
                                 temporalVdofs,
                                 heatFluxAtInitialTime);

    return heatFluxAtInitialTime;
}

Vector heat::SolutionHandler
:: getHeatFluxDataAtEndTime() const
{
    auto temporalMesh = m_temporalFeSpace->GetMesh();
    auto heatFluxData = getHeatFluxData();

    Array<int> temporalVdofs;
    const FiniteElement *temporalFe
            = m_temporalFeSpace->GetFE(temporalMesh->GetNE()-1);
    m_temporalFeSpace->GetElementVDofs(temporalMesh->GetNE()-1,
                                       temporalVdofs);

    Vector temporalShape(temporalVdofs.Size());
    IntegrationPoint ip;
    ip.Set1w(1.0, 1.0);
    temporalFe->CalcShape(ip, temporalShape);

    Vector heatFluxAtEndTime(m_spatialFeSpaceSizeForHeatFlux);
    buildSolutionAtSpecifiedTime(heatFluxData,
                                 temporalShape,
                                 temporalVdofs,
                                 heatFluxAtEndTime);

    return heatFluxAtEndTime;
}
