#pragma once

#include <filesystem>

#include "equation.h"

namespace enso {

// integrate implementation for the variational problem
DElNinoC1Solution var_integrate(
    DElNinoSetup& setup,
    int iters,
    DElNinoC1Solution const& initial,
    DElNinoC1Solution& result
);

std::pair<capd::DVector, capd::DMatrix> phiT(
    DElNinoSetup& setup,
    capd::DVector const& x,
    DElNinoSetup::TimePoint const& T
);

void draw(
    DElNinoSetup& setup,
    DElNinoSolution const& solution,
    int iters,
    std::filesystem::path const& output_dir,
    std::string const& name
);

// Use Newton's method to improve a periodic solution
DElNinoSolution newton_improve_periodic_solution(
    DElNinoSetup& setup,
    DElNinoSolution const& initial,
    DElNinoSetup::TimePoint T,
    double newton_value_precision,
    double newton_step_precision,
    double newton_maximum_iterations,
    std::ostream* output_ptr
);

// select C with first columns corresponding to most significant eigenvectors
// and the rest orthogonal to most significant left eigenvectors
capd::DMatrix most_significant_eigenvectors_left_orthogonal(
    DElNinoSetup& setup,
    DElNinoSolution const& X,
    DElNinoSetup::TimePoint const& T,
    int vectors_to_use,
    std::ostream* output_ptr
);

// select C with first columns corresponding to most significant eigenvectors
// and he rest orthogonal to most significant eigenvectors
capd::DMatrix most_significant_eigenvectors_orthogonal(
    DElNinoSetup& setup,
    DElNinoSolution const& X,
    DElNinoSetup::TimePoint const& T,
    int vectors_to_use,
    std::ostream* output_ptr
);

// select C with first columns corresponding to most significant eigenvectors
// and the rest from standard basis
capd::DMatrix most_significant_eigenvectors_identity(
    DElNinoSetup& setup,
    DElNinoSolution const& X,
    DElNinoSetup::TimePoint const& T,
    int vectors_to_use,
    std::ostream* output_ptr
);

capd::DMatrix identity(DElNinoSetup& setup, DElNinoSolution const&);

} // namespace enso
