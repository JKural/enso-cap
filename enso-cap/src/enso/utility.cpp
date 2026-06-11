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
} // namespace enso
