#include "enso-cap/enso/nonrigorous.h"

#include <complex>

#include "capd/alglib/lib.h"
#include "enso-cap/enso/utility.h"

namespace enso {

DElNinoC1Solution var_integrate(
    DElNinoSetup& setup,
    int iters,
    DElNinoC1Solution const& initial,
    DElNinoC1Solution& result
) {
    // create solver for the variational equation
    auto solver = setup.makeC1Solver(setup.makeC1Equation());

    // a throwaway vector of jacobians, they are stored inside the solution
    // anyways
    DElNinoSetup::Jacobians Js;

    // solution over the entire integration time
    DElNinoC1Solution solution = initial;

    // integrate for iters steps
    for (int i = 0; i < iters; ++i) {
        solver(solution, Js);
    }

    // solution over [currentTime - tau, currentTime]
    DElNinoC1Solution tail =
        solution.subcurve(solution.currentTime() - setup.tau());

    // copy the tail into result with correct taylor polynomial degrees
    capd::ddeshelper::copyReduce(tail, result);

    return solution;
}

std::pair<capd::DVector, capd::DMatrix> phiT(
    DElNinoSetup& setup,
    capd::DVector const& x,
    DElNinoSetup::TimePoint const& T
) {
    assert(setup.M() == x.dimension());
    // create solver for the variational equation
    auto solver = setup.makeC1Solver(setup.makeC1Equation());

    // placeholder section, used only to create a poincare map
    auto section = DElNinoSetup::C1JetSection(1, 0);

    // placeholder poincare map, used only to initialize the variational matrix to
    // Id
    auto poincare_map = setup.makeC1PoincareMap(solver, section);

    // prepare a solution for integration
    auto X = setup.makeC1Solution(x);

    // initialize derivative matrices in the solution
    poincare_map.setInitialV(X);

    // create a result solution
    auto Y = setup.makeC1Solution();

    // matrix to store the derivative
    capd::DMatrix D(setup.M(), setup.M());

    // integrate the solution
    var_integrate(setup, int(T), X, Y);

    // extract the derivative
    solver.extractVariationalMatrix(Y, D);

    return std::pair {Y.get_x(), D};
}

void draw(
    DElNinoSetup& setup,
    DElNinoSolution const& solution,
    int iters,
    std::filesystem::path const& output_dir,
    std::string const& name
) {
    auto extended_solution = setup.integrate(iters, solution);
    setup.drawSolution(output_dir, name, extended_solution);
    setup.drawDelayMap(output_dir, name + "_delay_map", extended_solution);
}

DElNinoSolution newton_improve_periodic_solution(
    DElNinoSetup& setup,
    DElNinoSolution const& initial,
    DElNinoSetup::TimePoint T,
    double newton_value_precision,
    double newton_step_precision,
    double newton_maximum_iterations,
    std::ostream* output_ptr
) {
    // identity matrix
    capd::DMatrix const id = capd::DMatrix::Identity(setup.M());

    // x - phiT(x), used with Netwon's method to detect a fixed point
    auto f = [&setup, &id, &T](auto const& x) {
        auto [phiTx, DphiTx] = phiT(setup, x, T);
        return std::pair {x - phiTx, id - DphiTx};
    };

    // a function to measure how "big" x is
    auto measure = [](auto&& x) { return max_norm(x); };

    // essentialy an alias, necessary to fix ambiguity
    auto solve = [](auto&& A, auto&& y) {
        return capd::matrixAlgorithms::gauss(
            std::forward<decltype(A)>(A),
            std::forward<decltype(y)>(y)
        );
    };

    auto x = initial.get_x();

    auto x_improved = newton(
        f,
        x,
        solve,
        measure,
        newton_value_precision,
        newton_step_precision,
        newton_maximum_iterations,
        output_ptr
    );

    return setup.makeSolution(x_improved);
}

capd::DMatrix most_significant_eigenvectors_left_orthogonal(
    DElNinoSetup& setup,
    DElNinoSolution const& X,
    DElNinoSetup::TimePoint const& T,
    int vectors_to_use,
    std::ostream* output_ptr
) {
    // D is equal to d/dx phi(T, x)
    auto [y, D] = phiT(setup, X.get_x(), T);
    auto D_T = D;
    D_T.transpose();

    capd::DVector re_lambda(setup.M());
    capd::DVector im_lambda(setup.M());
    capd::DMatrix re_eigenvectors(setup.M(), setup.M());
    capd::DMatrix im_eigenvectors(setup.M(), setup.M());

    capd::DVector re_T_lambda(setup.M());
    capd::DVector im_T_lambda(setup.M());
    capd::DMatrix re_T_eigenvectors(setup.M(), setup.M());
    capd::DMatrix im_T_eigenvectors(setup.M(), setup.M());
    capd::computeEigenvaluesAndEigenvectors(
        D,
        re_lambda,
        im_lambda,
        re_eigenvectors,
        im_eigenvectors
    );
    capd::computeEigenvaluesAndEigenvectors(
        D_T,
        re_T_lambda,
        im_T_lambda,
        re_T_eigenvectors,
        im_T_eigenvectors
    );

    if (output_ptr) {
        using namespace std::complex_literals;
        *output_ptr << "Most significant eigenvalues:\n";
        for (int i = 0; i < vectors_to_use; ++i) {
            auto complex_T_lambda_i = re_T_lambda[i] + im_T_lambda[i] * 1i;
            *output_ptr << i + 1 << ". ";
            *output_ptr << complex_T_lambda_i;
            *output_ptr << ", abs: ";
            *output_ptr << std::abs(complex_T_lambda_i);
            *output_ptr << "\n";
        }
        auto complex_T_lambda_vtu =
            re_T_lambda[vectors_to_use] + im_T_lambda[vectors_to_use] * 1i;
        *output_ptr << "First unused eigenvalue: ";
        *output_ptr << complex_T_lambda_vtu;
        *output_ptr << ", abs: ";
        *output_ptr << std::abs(complex_T_lambda_vtu);
        *output_ptr << "\n\n";
    };

    capd::DMatrix Q(setup.M(), setup.M());
    capd::DMatrix R(setup.M(), setup.M());
    capd::DMatrix C(setup.M(), setup.M());
    capd::DMatrix QR_input = capd::DMatrix::Identity(setup.M());
    for (int i = 0; i < vectors_to_use; ++i) {
        QR_input.column(i) = re_T_eigenvectors.column(i);
    }

    try {
        capd::matrixAlgorithms::QR_decompose(QR_input, Q, R);
        for (int i = 0; i < vectors_to_use; ++i) {
            re_eigenvectors.column(i).normalize();
            C.column(i) = re_eigenvectors.column(i);
        }
        for (std::size_t i = vectors_to_use; i < setup.M(); ++i) {
            C.column(i) = Q.column(i);
        }
    } catch (...) {
        std::stringstream ss;
        ss << "ERROR: Failed to orthogonalize the matrix, there might ";
        ss << "exist significant imaginary eigenvalues";
        throw std::logic_error(ss.str());
    }
    return C;
}

capd::DMatrix most_significant_eigenvectors_orthogonal(
    DElNinoSetup& setup,
    DElNinoSolution const& X,
    DElNinoSetup::TimePoint const& T,
    int vectors_to_use,
    std::ostream* output_ptr
) {
    // D is equal to d/dx phi(T, x)
    auto [y, D] = phiT(setup, X.get_x(), T);
    // we calculate the eigenvalue decomposition of D in order to find the most
    // significant eigenvalues and their respective eigenvectors
    capd::DVector re_lambda(setup.M());
    capd::DVector im_lambda(setup.M());
    capd::DMatrix re_eigenvectors(setup.M(), setup.M());
    capd::DMatrix im_eigenvectors(setup.M(), setup.M());
    capd::computeEigenvaluesAndEigenvectors(
        D,
        re_lambda,
        im_lambda,
        re_eigenvectors,
        im_eigenvectors
    );

    if (output_ptr) {
        using namespace std::complex_literals;
        *output_ptr << "Most significant eigenvalues:\n";
        for (int i = 0; i < vectors_to_use; ++i) {
            auto complex_lambda_i = re_lambda[i] + im_lambda[i] * 1i;
            *output_ptr << i + 1 << ". ";
            *output_ptr << complex_lambda_i;
            *output_ptr << ", abs: ";
            *output_ptr << std::abs(complex_lambda_i);
            *output_ptr << "\n";
        }
        auto complex_lambda_vtu =
            re_lambda[vectors_to_use] + im_lambda[vectors_to_use] * 1i;
        *output_ptr << "First unused eigenvalue: ";
        *output_ptr << complex_lambda_vtu;
        *output_ptr << ", abs: ";
        *output_ptr << std::abs(complex_lambda_vtu);
        *output_ptr << "\n\n";
    };

    // We define C to be: first vectors_to_use eigenvectors from D
    // and then fill the rest with some orthonormal vectors
    capd::DMatrix Q(setup.M(), setup.M());
    capd::DMatrix R(setup.M(), setup.M());
    capd::DMatrix C(setup.M(), setup.M());
    capd::DMatrix QR_input = capd::DMatrix::Identity(setup.M());
    for (int i = 0; i < vectors_to_use; ++i) {
        QR_input.column(i) = re_eigenvectors.column(i);
    }

    try {
        capd::matrixAlgorithms::QR_decompose(QR_input, Q, R);
        for (int i = 0; i < vectors_to_use; ++i) {
            re_eigenvectors.column(i).normalize();
            C.column(i) = re_eigenvectors.column(i);
        }
        for (std::size_t i = vectors_to_use; i < setup.M(); ++i) {
            C.column(i) = Q.column(i);
        }
    } catch (...) {
        std::stringstream ss;
        ss << "ERROR: Failed to orthogonalize the matrix, there might ";
        ss << "exist significant imaginary eigenvalues";
        throw std::logic_error(ss.str());
    }
    return C;
}

capd::DMatrix most_significant_eigenvectors_identity(
    DElNinoSetup& setup,
    DElNinoSolution const& X,
    DElNinoSetup::TimePoint const& T,
    int vectors_to_use,
    std::ostream* output_ptr
) {
    // D is equal to d/dx phi(T, x)
    auto [y, D] = phiT(setup, X.get_x(), T);
    // we calculate the eigenvalue decomposition of D in order to find the most
    // significant eigenvalues and their respective eigenvectors
    capd::DVector re_lambda(setup.M());
    capd::DVector im_lambda(setup.M());
    capd::DMatrix re_eigenvectors(setup.M(), setup.M());
    capd::DMatrix im_eigenvectors(setup.M(), setup.M());
    capd::computeEigenvaluesAndEigenvectors(
        D,
        re_lambda,
        im_lambda,
        re_eigenvectors,
        im_eigenvectors
    );

    if (output_ptr) {
        using namespace std::complex_literals;
        *output_ptr << "Most significant eigenvalues:\n";
        for (int i = 0; i < vectors_to_use; ++i) {
            auto complex_lambda_i = re_lambda[i] + im_lambda[i] * 1i;
            *output_ptr << i + 1 << ". ";
            *output_ptr << complex_lambda_i;
            *output_ptr << ", abs: ";
            *output_ptr << std::abs(complex_lambda_i);
            *output_ptr << "\n";
        }
        auto complex_lambda_vtu =
            re_lambda[vectors_to_use] + im_lambda[vectors_to_use] * 1i;
        *output_ptr << "First unused eigenvalue: ";
        *output_ptr << complex_lambda_vtu;
        *output_ptr << ", abs: ";
        *output_ptr << std::abs(complex_lambda_vtu);
        *output_ptr << "\n\n";
    };

    // We define C to be: first vectors_to_use eigenvectors from D
    // and then fill the rest with vectors from canonical basis
    capd::DMatrix C(setup.M(), setup.M());
    capd::DMatrix Id = capd::DMatrix::Identity(setup.M());
    for (int i = 0; i < vectors_to_use; ++i) {
        re_eigenvectors.column(i).normalize();
        C.column(i) = re_eigenvectors.column(i);
    }
    for (std::size_t i = vectors_to_use; i < setup.M(); ++i) {
        C.column(i) = Id.column(i);
    }
    return C;
}

capd::DMatrix identity(DElNinoSetup& setup, DElNinoSolution const&) {
    return capd::DMatrix::Identity(setup.M());
}

} // namespace enso
