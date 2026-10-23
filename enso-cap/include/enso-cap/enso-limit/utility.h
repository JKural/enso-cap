#pragma once

#include "equation.h"

namespace enso {

std::pair<DElNinoLimitSetup, IElNinoLimitSetup> create_setups_limit(
    capd::interval alpha,
    capd::interval beta,
    capd::interval tau,
    int p,
    int n,
    int reqSteps = 0,
    int maxSteps = 0,
    unsigned int max_order = 99
);

std::string
tabular_print(DElNinoLimitSetup& setup, DElNinoLimitSolution const& x);

std::string
tabular_print(IElNinoLimitSetup& setup, IElNinoLimitSolution const& x);

} // namespace enso
