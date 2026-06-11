#include <enso-cap/enso-cap.h>

using namespace enso;

int main() {
    double alpha = 1.0;
    double beta = 3.0;
    double kappa = 11.0;
    double tau = 1.2;
    int p = 100;
    int n = 10;
    double T0 = 1000;
    std::string output_directory = ".";
    std::string output_name = "plot";
    int integration_iters = static_cast<int>(T0 / tau * p);
    capd::DMap x0map(
        "par:lambda;"
        "var:x;"
        "fun:lambda*x;",
        n + 1
    );
    DElNinoSetup
        setup(DElNinoParamsVector {alpha, beta, kappa, tau}, p, n, 0, 0, n);
    double lambda_min = 0;
    double lambda_max = 1;
    for (auto i = 0; i < 100; ++i) {
        double lambda = (lambda_min + lambda_max) / 2;
        x0map.setParameter("lambda", lambda);
        DElNinoSolution
            x0(setup.grid().point(-p), setup.grid().point(0), n, x0map);
        auto x1 = setup.makeSolution();
        setup.integrate(integration_iters, x0, x1);
        x1 = setup.integrate(5 / tau * p, x1);
        double max = 1.4;
        bool skip = false;
        for (int j = 0; j < p; ++j) {
            if (abs(x1.value(setup.grid().point(-j))[0]) > max) {
                lambda_max = lambda;
                skip = true;
                break;
            }
        }
        if (skip) {
            continue;
        }
        lambda_min = lambda;
    }
    auto lambda = (lambda_max + lambda_min) / 2;
    std::cout << "lambda = " << lambda << "\n";

    auto eps = std::numeric_limits<double>::epsilon();
    x0map.setParameter("lambda", lambda - eps);
    DElNinoSolution
        x0_short(setup.grid().point(-p), setup.grid().point(0), n, x0map);
    draw(setup, x0_short, 1000 / tau * p, ".", "fixed-point-short");
    x0map.setParameter("lambda", lambda + eps);
    DElNinoSolution
        x0_long(setup.grid().point(-p), setup.grid().point(0), n, x0map);
    draw(setup, x0_long, 1000 / tau * p, ".", "fixed-point-long");
}
