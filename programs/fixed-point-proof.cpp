#include <enso-cap/enso-cap.h>

struct Config {
    capd::interval a = 1.0;
    capd::interval b = 3.0;
    capd::interval kappa = 0.34;
    capd::interval tau = capd::interval(12) / capd::interval(10);
    double x0 = 0.0;
    int p = 60;
    int n = 2;
    capd::interval T0 = 100.0;
    capd::interval T = 2.0;
    int newton_iters = 10;
    double newton_error = 1e-10;
    std::string output_directory = ".";
    std::string output_name = "fixed-point-proof";
};

std::string print(Config const& cfg) {
    std::stringstream sstream;
    sstream << "a                = " << cfg.a << "\n";
    sstream << "b                = " << cfg.b << "\n";
    sstream << "kappa            = " << cfg.kappa << "\n";
    sstream << "tau              = " << cfg.tau << "\n";
    sstream << "x0               = " << cfg.x0 << "\n";
    sstream << "p                = " << cfg.p << "\n";
    sstream << "n                = " << cfg.n << "\n";
    sstream << "T0               = " << cfg.T0 << "\n";
    sstream << "T                = " << cfg.T << "\n";
    sstream << "newton_iters     = " << cfg.newton_iters << "\n";
    sstream << "newton_error     = " << cfg.newton_error << "\n";
    sstream << "output_directory = " << cfg.output_directory << "\n";
    sstream << "output_name      = " << cfg.output_name << "\n";
    return sstream.str();
}

capd::ddeshelper::ArgumentParser
parse_commandline(char argc, char** argv, Config& cfg) {
    capd::ddeshelper::ArgumentParser args(argc, argv);
    args.parse("--a=", cfg.a, "Parameter a");
    args.parse("--b=", cfg.b, "Parameter b");
    args.parse("--kappa=", cfg.kappa, "Parameter kappa");
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
    args.parse("--output_directory=", cfg.output_directory, "Output directory");
    args.parse("--output_name=", cfg.output_name, "Output name");
    return args;
}

using namespace enso;

int main(int argc, char** argv) {
    Config cfg;
    auto args = parse_commandline(argc, argv, cfg);
    if (args.isHelpRequested()) {
        std::cout << args.getHelp() << "\n";
        return 0;
    }
    int init_integration_iters =
        static_cast<int>(midpoint(cfg.T0 / cfg.tau * cfg.p));
    auto [setup, isetup] =
        create_setups(cfg.a, cfg.b, cfg.kappa, cfg.tau, cfg.p, cfg.n, 0, 0, 99);
    DElNinoSolution x0(
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
        setup.t(midpoint(cfg.T / cfg.tau * cfg.p)),
        cfg.newton_error,
        cfg.newton_error,
        cfg.newton_iters,
        &std::cout
    );
    draw(
        setup,
        Nx,
        midpoint(cfg.T / cfg.tau * cfg.p),
        cfg.output_directory,
        cfg.output_name + "-candidate"
    );

    auto ix0 = cast_to_interval(
        isetup,
        setup,
        Nx,
        capd::interval(-1, 1) * 3e-3,
        capd::interval(-1, 1) * 1e-0,
        most_significant_eigenvectors_left_orthogonal,
        setup.grid().point(midpoint(cfg.T / cfg.tau * cfg.p)),
        1,
        &std::cout
    );
    auto C_inv = capd::matrixAlgorithms::gaussInverseMatrix(ix0.get_C());
    for (int i = 0; i < 5; ++i) {
        inflate(
            isetup,
            ix0,
            isetup.grid().point(midpoint(cfg.T / cfg.tau * cfg.p)),
            C_inv,
            1.05,
            &std::cout
        );
    }

    auto [ix, r1] = split_integrate(
        isetup,
        ix0,
        isetup.grid().point(midpoint(cfg.T / cfg.tau * cfg.p)),
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

    draw(
        isetup,
        ix0,
        1,
        cfg.output_directory,
        cfg.output_name + "-starting-set"
    );
    draw(isetup, ix, 1, cfg.output_directory, cfg.output_name + "-result-set");
}
