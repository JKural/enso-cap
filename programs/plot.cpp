#include <enso-cap/enso-cap.h>

using namespace enso;

struct Config {
    double a = 1.0;
    double b = 3.0;
    double kappa = 11.0;
    double tau = 1.2;
    double x0 = 0.5;
    int p = 100;
    int n = 10;
    double T0 = 100.0;
    double T1 = 110.0;
    std::string output_directory = ".";
    std::string output_name = "plot";
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
    sstream << "T1               = " << cfg.T1 << "\n";
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
    args.parse("--T1=", cfg.T1, "Final integration time");
    args.parse("--output_directory=", cfg.output_directory, "Output directory");
    args.parse("--output_name=", cfg.output_name, "Output name");
    return args;
}

int main(int argc, char** argv) {
    Config cfg;
    auto args = parse_commandline(argc, argv, cfg);
    if (args.isHelpRequested()) {
        std::cout << args.getHelp() << "\n";
        return 0;
    }
    std::cout << print(cfg);
    int init_integration_iters = static_cast<int>(cfg.T0 / cfg.tau * cfg.p);
    int drawing_integration_iters =
        static_cast<int>((cfg.T1 - cfg.T0) / cfg.tau * cfg.p);
    DElNinoSetup setup(
        DElNinoParamsVector {cfg.a, cfg.b, cfg.kappa, cfg.tau},
        cfg.p,
        cfg.n,
        0,
        0,
        cfg.n
    );
    DElNinoSolution x0(
        setup.grid().point(-cfg.p),
        setup.grid().point(0),
        cfg.n,
        capd::DVector {cfg.x0}
    );
    auto x1 = setup.makeSolution();
    setup.integrate(init_integration_iters, x0, x1);
    draw(
        setup,
        x1,
        drawing_integration_iters,
        cfg.output_directory,
        cfg.output_name
    );
}
