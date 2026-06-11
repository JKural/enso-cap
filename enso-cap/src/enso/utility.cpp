#include "enso-cap/enso/utility.h"

namespace enso {

std::pair<DElNinoSetup, IElNinoSetup> create_setups(
    capd::interval alpha,
    capd::interval beta,
    capd::interval kappa,
    capd::interval tau,
    int p,
    int n,
    int req_steps,
    int max_steps,
    unsigned int max_order
) {
    return {
        DElNinoSetup(
            DElNinoSetup::ParamsVector {
                midpoint(alpha),
                midpoint(beta),
                midpoint(kappa),
                midpoint(tau),
            },
            p,
            n,
            req_steps,
            max_steps,
            max_order
        ),
        IElNinoSetup(
            p,
            n,
            IElNinoSetup::ParamsVector {
                alpha,
                beta,
                kappa,
                tau,
            },
            req_steps,
            max_steps,
            max_order
        )
    };
}

std::string tabular_print(DElNinoSetup& setup, DElNinoSolution const& X) {
    std::stringstream sstream;
    sstream << std::setprecision(17);
    auto const& x = X.get_x();
    auto n = setup.n();
    sstream << "t";
    for (unsigned i = 0; i <= n; ++i) {
        sstream << ", f^(" << i << ")";
    }
    sstream << "\n";
    assert((x.dimension() - 1) % (n + 1) == 0);
    for (unsigned i = 1; i < x.dimension() / (n + 1); ++i) {
        sstream << static_cast<double>(
            X.getCurrentTime() - setup.grid().point(i)
        );
        for (unsigned j = 0; j <= n; ++j) {
            sstream << ", " << x[1 + (n + 1) * i + j];
        }
        sstream << "\n";
    }
    return sstream.str();
}

std::string tabular_print(IElNinoSolution const& x);

} // namespace enso
