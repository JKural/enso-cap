#pragma once

#include "equation.h"

namespace enso {

template<class T>
typename T::ScalarType max_norm(T const& v) {
    return capd::abs(
        *std::max_element(v.begin(), v.end(), [](auto&& x, auto&& y) {
            return capd::abs(x) < capd::abs(y);
        })
    );
}

// Generic implementation of Newton's method
template<class F, class S, class M, class T>
T newton(
    F f,
    T x,
    S solve,
    M measure,
    double val_eps,
    double step_eps,
    int max_iters,
    std::ostream* output_ptr
) {
    if (output_ptr) {
        *output_ptr << "Starting Newton's method\n";
        *output_ptr << "val_eps =   " << val_eps << "\n";
        *output_ptr << "step_eps =  " << step_eps << "\n";
        *output_ptr << "max_iters = " << max_iters << "\n";
        *output_ptr << "\n";
    }

    int i = 0;
    double val_diff = std::numeric_limits<double>::infinity();
    double step_diff = std::numeric_limits<double>::infinity();

    for (; i < max_iters; ++i) {
        auto [fx, dfx] = f(x);
        auto nx = solve(dfx, dfx * x - fx);
        val_diff = measure(fx);
        step_diff = measure(nx - x);
        x = nx;
        if (output_ptr) {
            *output_ptr << "Step " << i + 1 << "/" << max_iters << ":\n";
            *output_ptr << "val_diff =  " << val_diff << "\n";
            *output_ptr << "step_diff = " << step_diff << "\n";
            *output_ptr << "\n";
        }
        if (val_diff < val_eps || step_diff < step_eps) {
            break;
        }
    }

    if (output_ptr) {
        *output_ptr << "Finished Newton's method\n";
        *output_ptr << "Final val_diff =  " << val_diff << ", ";
        *output_ptr << "val_eps =  " << val_eps << "\n";
        *output_ptr << "Final step_diff = " << step_diff << ", ";
        *output_ptr << "step_eps = " << step_eps << "\n";
        if (val_diff >= val_eps && step_diff >= step_eps) {
            *output_ptr << "WARNING: Failed to reach specified ";
            *output_ptr << "accuracy on function value or ";
            *output_ptr << "argument step";
        }
        *output_ptr << "\n";
    }

    return x;
}

template<class Interval>
double midpoint(Interval i) {
    return (i.leftBound() + i.rightBound()) / 2.0;
}

std::pair<DElNinoSetup, IElNinoSetup> create_setups(
    capd::interval alpha,
    capd::interval beta,
    capd::interval kappa,
    capd::interval tau,
    int p,
    int n,
    int reqSteps = 0,
    int maxSteps = 0,
    unsigned int max_order = 99
);

std::string tabular_print(DElNinoSetup& setup, DElNinoSolution const& x);

std::string tabular_print(IElNinoSetup& setup, IElNinoSolution const& x);

} // namespace enso
