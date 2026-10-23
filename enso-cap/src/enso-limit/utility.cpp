#include "enso-cap/enso-limit/utility.h"

#include "enso-cap/enso/utility.h"

namespace enso {

std::pair<DElNinoLimitSetup, IElNinoLimitSetup> create_setups_limit(
    capd::interval alpha,
    capd::interval beta,
    capd::interval tau,
    int p,
    int n,
    int req_steps,
    int max_steps,
    unsigned int max_order
) {
    return {
        DElNinoLimitSetup(
            DElNinoLimitSetup::ParamsVector {
                midpoint(alpha),
                midpoint(beta),
                midpoint(tau)
            },
            p,
            n,
            req_steps,
            max_steps,
            max_order
        ),
        IElNinoLimitSetup(
            p,
            n,
            IElNinoLimitSetup::ParamsVector {
                alpha,
                beta,
                tau,
            },
            req_steps,
            max_steps,
            max_order
        )
    };
}

std::string
tabular_print(DElNinoLimitSetup& setup, DElNinoLimitSolution const& X) {
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

std::string tabular_print(IElNinoLimitSolution const& x);

} // namespace enso
