#include "enso-cap/enso/rigorous.h"

#include <numeric>

#include "enso-cap/enso/utility.h"

namespace enso {

void draw(
    IElNinoSetup& setup,
    IElNinoSolution const& solution,
    int iters,
    std::filesystem::path const& output_dir,
    std::string const& name
) {
    auto solution_copy = solution;
    auto extended_solution = setup.timemap(solution_copy, iters);
    setup.drawSolution(output_dir, name, extended_solution);
    // rigorous helper doesnt implement drawing delay maps
    // setup.drawDelayMap(output_dir, name + "_delay_map", extended_solution);
}

void inflate(
    IElNinoSetup& setup,
    IElNinoSolution& X,
    IElNinoSetup::TimePoint const& T,
    capd::IMatrix const& C_inv,
    capd::interval inflation_coefficient,
    std::ostream* output_ptr
) {
    auto X_copy = X;
    auto r0 = X.get_r0();
    auto xi = X.get_Xi();
    double r0_norm_before = max_norm(r0).rightBound();
    double xi_norm_before = max_norm(xi).rightBound();

    // push X by one period
    auto Y = setup.timemap(X_copy, int(T));

    // choose next r0 and xi such that they contain previous and new values
    r0 = capd::vectalg::intervalHull(
        r0,
        C_inv * Y.get_x() - C_inv * X.get_x() + C_inv * Y.get_C() * Y.get_r0()
            + C_inv * Y.get_B() * Y.get_r()
    );
    xi = capd::vectalg::intervalHull(xi, Y.get_Xi());

    // symmetrize r0
    r0 *= capd::interval {-1, 1};

    // further expand r0 and xi
    r0 *= inflation_coefficient;

    // xi is not centered, so we need additional work
    auto center = xi;
    capd::vectalg::mid(xi, center);
    xi -= center;
    xi *= inflation_coefficient;
    xi += center;

    if (output_ptr) {
        double r0_norm_after = max_norm(r0).rightBound();
        double xi_norm_after = max_norm(xi).rightBound();

        *output_ptr << "r0 norm before inflation: " << r0_norm_before << "\n";
        *output_ptr << "r0 norm after  inflation: " << r0_norm_after << "\n";
        *output_ptr << "xi norm before inflation: " << xi_norm_before << "\n";
        *output_ptr << "xi norm after  inflation: " << xi_norm_after << "\n\n";
    }

    X.set_r0(r0);
    X.set_Xi(xi);
}

std::pair<IElNinoSolution, capd::IVector> split_integrate(
    IElNinoSetup& setup,
    IElNinoSolution const& X,
    IElNinoSetup::TimePoint const& T,
    int n_of_splits,
    capd::IMatrix const& C_inv
) {
    auto delta_r = (X.get_r0()[0].right() - X.get_r0()[0].left()) / n_of_splits;
    auto delta_r_interval = intervalHull(0, delta_r);
    std::vector<int> splits(n_of_splits);
    // vector of split solutions and their r0 vectors in original variables
    std::vector<std::pair<IElNinoSolution, capd::IVector>> results(
        n_of_splits,
        std::pair {IElNinoSolution(setup.grid().point(0)), capd::IVector()}
    );
    std::iota(splits.begin(), splits.end(), 0);
    auto integrate = [&](int i) -> std::pair<IElNinoSolution, capd::IVector> {
        auto r0_split = X.get_r0();
        r0_split[0] = r0_split[0].left() + i * delta_r + delta_r_interval;
        auto x_split = X.get_x();
        auto center = capd::vectalg::midVector(r0_split);
        r0_split -= center;
        x_split += X.get_C() * center;
        // capd::interval center = r0_split[0].mid();
        // r0_split[0] -= center;
        // x_split += X.get_C().column(1) * center;
        IElNinoSolution X_split = X;
        X_split.set_r0(r0_split);
        X_split.set_x(x_split);
        auto X_split_copy = X_split;
        auto Y_split = setup.timemap(X_split_copy, static_cast<int>(T));
        auto r1_split = C_inv * Y_split.get_x() - C_inv * X_split.get_x()
            + C_inv * Y_split.get_C() * Y_split.get_r0()
            + C_inv * Y_split.get_B() * Y_split.get_r()
            /*- C_inv * X_split.get_B() * X_split.get_r()*/;
        return {Y_split, r1_split};
    };
    auto merge = [&](
                     std::pair<IElNinoSolution, capd::IVector> const& lhs,
                     std::pair<IElNinoSolution, capd::IVector> const& rhs
                 ) -> std::pair<IElNinoSolution, capd::IVector> {
        auto result_solution = lhs.first;
        result_solution.set_x(
            intervalHull(lhs.first.get_x(), rhs.first.get_x())
        );
        result_solution.set_C(
            intervalHull(lhs.first.get_C(), rhs.first.get_C())
        );
        result_solution.set_r0(
            intervalHull(lhs.first.get_r0(), rhs.first.get_r0())
        );
        result_solution.set_Xi(
            intervalHull(lhs.first.get_Xi(), rhs.first.get_Xi())
        );
        result_solution.set_B(
            intervalHull(lhs.first.get_B(), rhs.first.get_B())
        );
        result_solution.set_r(
            intervalHull(lhs.first.get_r(), rhs.first.get_r())
        );
        auto result_r1 = intervalHull(lhs.second, rhs.second);
        return std::pair {result_solution, result_r1};
    };
    std::transform(splits.begin(), splits.end(), results.begin(), integrate);
    return std::reduce(results.begin() + 1, results.end(), results[0], merge);
}

} // namespace enso
