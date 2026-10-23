#include <enso-cap/enso-cap.h>

using namespace enso;

struct Config {
    double a = 1.0;
    double b = 3.0;
    double tau = 1.2;
    double x0 = 0.0;
    int p = 60;
    int n = 2;
    double T0 = 100.0;
    double T = 1.0;
    int newton_iters = 10;
    double newton_error = 1e-10;
};

std::string print(Config const& cfg) {
    std::stringstream sstream;
    sstream << "a                = " << cfg.a << "\n";
    sstream << "b                = " << cfg.b << "\n";
    sstream << "tau              = " << cfg.tau << "\n";
    sstream << "x0               = " << cfg.x0 << "\n";
    sstream << "p                = " << cfg.p << "\n";
    sstream << "n                = " << cfg.n << "\n";
    sstream << "T0               = " << cfg.T0 << "\n";
    sstream << "T                = " << cfg.T << "\n";
    sstream << "newton_iters     = " << cfg.newton_iters << "\n";
    sstream << "newton_error     = " << cfg.newton_error << "\n";
    return sstream.str();
}

capd::ddeshelper::ArgumentParser
parse_commandline(char argc, char** argv, Config& cfg) {
    capd::ddeshelper::ArgumentParser args(argc, argv);
    args.parse("--a=", cfg.a, "Parameter a");
    args.parse("--b=", cfg.b, "Parameter b");
    args.parse("--tau=", cfg.tau, "Delay");
    args.parse("--x0=", cfg.x0, "Value of the initial (constant) function");
    args.parse("--p=", cfg.p, "Number of subintervals per delay");
    args.parse("--n=", cfg.n, "Order of the method");
    args.parse("--T0=", cfg.T0, "Initial integration time");
    args.parse("--T=", cfg.T, "Period");
    args.parse(
        "--newton_iters=",
        cfg.newton_iters,
        "Maximum number of Newton's method iterations"
    );
    args.parse(
        "--newton_error=",
        cfg.newton_error,
        "Target error for the Netwon's method"
    );
    return args;
}

int main(int argc, char** argv) {
    Config cfg;
    auto args = parse_commandline(argc, argv, cfg);
    if (args.isHelpRequested()) {
        std::cout << args.getHelp() << "\n";
        return 0;
    }
    std::cout << print(cfg) << "\n";
    int init_integration_iters = static_cast<int>(cfg.T0 / cfg.tau * cfg.p);
    DElNinoLimitSetup setup(
        DElNinoLimitParamsVector {cfg.a, cfg.b, cfg.tau},
        cfg.p,
        cfg.n,
        0,
        0,
        cfg.n
    );
    DElNinoLimitSolution x0(
        setup.grid().point(-cfg.p),
        setup.grid().point(0),
        cfg.n,
        capd::DVector {cfg.x0}
    );
    auto x1 = setup.makeSolution();
    setup.integrate(init_integration_iters, x0, x1);
    auto Nx = newton_improve_periodic_solution(
        setup,
        x0,
        setup.t(cfg.T / cfg.tau * cfg.p),
        cfg.newton_error,
        cfg.newton_error,
        cfg.newton_iters,
        &std::cout
    );
    auto [y, D] = phiT(setup, Nx.get_x(), setup.t(cfg.T / cfg.tau * cfg.p));
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
    for (std::size_t i = 0; i < re_lambda.dimension(); ++i) {
        std::cout << i + 1 << ". re: " << re_lambda[i]
                  << ", im: " << im_lambda[i] << "\n";
    }
}
