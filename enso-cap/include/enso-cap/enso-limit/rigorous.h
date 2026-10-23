#pragma once

#include <filesystem>

#include "equation.h"

namespace enso {

void draw(
    IElNinoLimitSetup& setup,
    IElNinoLimitSolution const& solution,
    int iters,
    std::filesystem::path const& output_dir,
    std::string const& name
);

// Cast to rigorous solution with a given choice strategy for matrix C
template<class F, class... Args>
IElNinoLimitSolution cast_to_interval_limit(
    IElNinoLimitSetup& isetup,
    DElNinoLimitSetup& setup,
    DElNinoLimitSolution const& X,
    capd::interval r0,
    capd::interval xi,
    F C_choice_strategy,
    Args&&... args
) {
    /*RSolution rX(rsetup.grid(), X.get_x().dimension());*/
    IElNinoLimitSolution iX(
        isetup.grid().point(-isetup.p()),
        isetup.grid().point(0),
        isetup.n(),
        capd::IVector {0.0}
    );
    iX.set_x(static_cast<capd::IVector>(X.get_x()));
    capd::IMatrix C = static_cast<capd::IMatrix>(
        C_choice_strategy(setup, X, std::forward<Args>(args)...)
    );
    capd::IVector r0_full_dimensional(C.numberOfColumns());
    for (auto& x : r0_full_dimensional) {
        x = r0;
    }
    capd::IVector xi_full_dimensional = iX.get_Xi();
    for (auto& x : xi_full_dimensional) {
        x = xi;
    }

    iX.set_Cr0(C, r0_full_dimensional);
    iX.set_Xi(xi_full_dimensional);

    return iX;
}

// inflate r0 and Xi
void inflate(
    IElNinoLimitSetup& setup,
    IElNinoLimitSolution& X,
    IElNinoLimitSetup::TimePoint const& T,
    capd::IMatrix const& C_inv,
    capd::interval inflation_coefficient_r0,
    capd::interval inflation_coefficient_xi,
    std::ostream* output_ptr
);

std::pair<IElNinoLimitSolution, capd::IVector> split_integrate(
    IElNinoLimitSetup& setup,
    IElNinoLimitSolution const& X,
    IElNinoLimitSetup::TimePoint const& T,
    int n_of_splits,
    capd::IMatrix const& C_inv
);

} // namespace enso
