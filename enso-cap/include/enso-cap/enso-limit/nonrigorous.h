#pragma once

#include <filesystem>

#include "equation.h"

namespace enso {

// integrate implementation for the variational problem
DElNinoLimitC1Solution var_integrate(
    DElNinoLimitSetup& setup,
    int iters,
    DElNinoLimitC1Solution const& initial,
    DElNinoLimitC1Solution& result
);

std::pair<capd::DVector, capd::DMatrix> phiT(
    DElNinoLimitSetup& setup,
    capd::DVector const& x,
    DElNinoLimitSetup::TimePoint const& T
);

void draw(
    DElNinoLimitSetup& setup,
    DElNinoLimitSolution const& solution,
    int iters,
    std::filesystem::path const& output_dir,
    std::string const& name
);

// Use Newton's method to improve a periodic solution
DElNinoLimitSolution newton_improve_periodic_solution(
    DElNinoLimitSetup& setup,
    DElNinoLimitSolution const& initial,
    DElNinoLimitSetup::TimePoint T,
    double newton_value_precision,
    double newton_step_precision,
    double newton_maximum_iterations,
    std::ostream* output_ptr
);

// select C with first columns corresponding to most significant eigenvectors
// and the rest orthogonal to most significant left eigenvectors
capd::DMatrix most_significant_eigenvectors_left_orthogonal_limit(
    DElNinoLimitSetup& setup,
    DElNinoLimitSolution const& X,
    DElNinoLimitSetup::TimePoint const& T,
    int vectors_to_use,
    std::ostream* output_ptr
);

// select C with first columns corresponding to most significant eigenvectors
// and he rest orthogonal to most significant eigenvectors
capd::DMatrix most_significant_eigenvectors_orthogonal_limit(
    DElNinoLimitSetup& setup,
    DElNinoLimitSolution const& X,
    DElNinoLimitSetup::TimePoint const& T,
    int vectors_to_use,
    std::ostream* output_ptr
);

// select C with first columns corresponding to most significant eigenvectors
// and the rest from standard basis
capd::DMatrix most_significant_eigenvectors_identity_limit(
    DElNinoLimitSetup& setup,
    DElNinoLimitSolution const& X,
    DElNinoLimitSetup::TimePoint const& T,
    int vectors_to_use,
    std::ostream* output_ptr
);

capd::DMatrix
identity_limit(DElNinoLimitSetup& setup, DElNinoLimitSolution const&);

} // namespace enso
