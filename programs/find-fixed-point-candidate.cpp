#include <enso-cap/enso-cap.h>

using namespace enso;

int main() {
    double alpha = 1.0;
    double beta = 3.0;
    double kappa = 11.0;
    double tau = 1.2;
    int p = 100;
    int n = 10;
    double T0 = 100;
    double period = 5;
    std::string output_directory = ".";
    int init_integration_iters = static_cast<int>(T0 / tau * p);
    capd::DMap x0map(
        "var:x;"
        "fun:x;",
        n + 1
    );
    DElNinoSetup
        setup(DElNinoParamsVector {alpha, beta, kappa, tau}, p, n, 0, 0, n);
    DElNinoSolution x0(setup.grid().point(-p), setup.grid().point(0), n, x0map);
    auto x1 = setup.makeSolution();
    setup.integrate(init_integration_iters, x0, x1);
    draw(setup, x1, 1, ".", "pre-newton");
    auto Nx = newton_improve_periodic_solution(
        setup,
        x0,
        setup.t(period / tau * p),
        1e-15,
        1e-15,
        100,
        &std::cout
    );
    draw(setup, Nx, 1, ".", "post-newton");
}
