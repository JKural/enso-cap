#include <enso-cap/enso-cap.h>

using namespace enso;

int main() {
    capd::interval alpha(1.0);
    capd::interval beta(3.0);
    capd::interval kappa(0.34);
    capd::interval tau = capd::interval(12) / capd::interval(10);
    capd::interval period(2);
    int p = 80;
    int n = 2;
    double T0 = 100;
    std::string output_directory = ".";
    int init_integration_iters = static_cast<int>(midpoint(T0 / tau * p));
    capd::DMap x0map(
        "var:x;"
        "fun:x;",
        n + 1
    );
    auto [setup, isetup] =
        create_setups(alpha, beta, kappa, tau, p, n, 0, 0, 10);
    DElNinoSolution x0(setup.grid().point(-p), setup.grid().point(0), n, x0map);
    auto x1 = setup.makeSolution();
    setup.integrate(init_integration_iters, x0, x1);
    auto Nx = newton_improve_periodic_solution(
        setup,
        x0,
        setup.t(midpoint(period / tau * p)),
        1e-15,
        1e-15,
        10,
        &std::cout
    );
    draw(setup, Nx, midpoint(period / tau * p), ".", "post-newton");

    auto ix0 = cast_to_interval(
        isetup,
        setup,
        Nx,
        capd::interval(-1, 1) * 3e-3,
        capd::interval(-1, 1) * 1e-0,
        most_significant_eigenvectors_left_orthogonal,
        setup.grid().point(midpoint(period / tau * p)),
        1,
        &std::cout
    );
    auto C_inv = capd::matrixAlgorithms::gaussInverseMatrix(ix0.get_C());
    for (int i = 0; i < 5; ++i) {
        inflate(
            isetup,
            ix0,
            isetup.grid().point(midpoint(period / tau * p)),
            C_inv,
            1.05,
            &std::cout
        );
    }

    auto [ix, r1] = split_integrate(
        isetup,
        ix0,
        isetup.grid().point(midpoint(period / tau * p)),
        8,
        C_inv
    );

    std::cout << std::setprecision(6);
    for (unsigned int i = 0; i < ix.get_r0().dimension(); ++i) {
        if (subsetInterior(r1[i], ix0.get_r0()[i])) {
            std::cout << "[X] ";
        } else {
            std::cout << "[ ] ";
        }
        std::cout << "r0[" << i << "]: " << r1[i] << " " << ix0.get_r0()[i]
                  << "\n";
    }

    for (unsigned int i = 0; i < ix.get_Xi().dimension(); ++i) {
        if (subsetInterior(ix.get_Xi()[i], ix0.get_Xi()[i])) {
            std::cout << "[X] ";
        } else {
            std::cout << "[ ] ";
        }
        std::cout << "Xi[" << i << "]: " << ix.get_Xi()[i] << " "
                  << ix0.get_Xi()[i] << "\n";
    }

    if (subsetInterior(r1, ix0.get_r0())
        && subsetInterior(ix.get_Xi(), ix0.get_Xi())) {
        std::cout << "Proof successfull\n";
    }

    draw(isetup, ix0, 1, ".", "ix0");
    draw(isetup, ix, 1, ".", "ix");
}
