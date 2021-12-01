#include "spatial_error_evaluator.hpp"

using namespace mfem;


double sparseHeat::SpatialErrorOfGradientOfTemperature
:: eval(std::shared_ptr<GridFunction> &,
        mfem::Array<double> &,
        mfem::Array<std::shared_ptr<GridFunction> > &) const
{
    return 0;
}
