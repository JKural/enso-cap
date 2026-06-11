#pragma once

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wreorder"
#pragma GCC diagnostic ignored "-Wsign-compare"
#include <capd/capdlib.h>
#include <capd/ddes/ddeslib.h>
#include <capd/ddeshelper/ddeshelperlib.h>

#include <capd/ddeshelper/DDEHelperNonrigorous.hpp>
#include <capd/ddeshelper/DDEHelperRigorous.hpp>
#pragma GCC diagnostic pop

namespace enso {

template<typename T>
inline T get_pi() {
    return T::pi();
}

template<>
inline double get_pi() {
    return 3.14159265358979323846;
}

// fadbad doesn't have tanh implemented - we use this formula instead
// could introduce errors
template<typename T>
inline T tanh(T const& x) {
    return (exp(x) - exp(-x)) / (exp(x) + exp(-x));
}

template<class Scalar, class Param = Scalar>
class ElNinoEquation {
public:
    using ScalarType = Scalar;
    using size_type = unsigned int;
    typedef ScalarType RealType;
    using ParamType = Param;
    using VectorType = capd::vectalg::Vector<ScalarType, 0>;
    using ParamsVectorType = capd::vectalg::Vector<ParamType, 0>;

    ElNinoEquation() : alpha(1), beta(1), kappa(1) {}

    ElNinoEquation(ParamType alpha, ParamType beta, ParamType kappa) :
        alpha(alpha),
        beta(beta),
        kappa(kappa) {}

    ElNinoEquation(capd::vectalg::Vector<ParamType, 0> const& params) :
        alpha(params[0]),
        beta(params[1]),
        kappa(params[2]) {}

    static size_type imageDimension() {
        return 1;
    };

    static size_type dimension() {
        return 2;
    }

    static size_type getParamsCount() {
        return 3;
    }

    template<typename Real, typename InVector, typename OutVector>
    void operator()(Real const& t, InVector const x, OutVector& fx) const {
        // a hack to make depreciation warning about implicit copy construction
        // go away
        // instead of 'auto Ttau = x[1]' we create an empty object and then copy
        // with operator=
        typename InVector::ScalarType Ttau;
        Ttau = x[1];
        auto const pi = get_pi<ScalarType>();
        fx[0] = -alpha * tanh(kappa * Ttau) + beta * cos(2 * pi * t);
    }

    static std::string show() {
        return "El Nino equation as autonomous system: $T'(t) = \\beta \\cos(2 \\pi s) - \\alpha \\tanh( \\kappa T(t-\\tau))$.";
    }

protected:
    ParamType alpha;
    ParamType beta;
    ParamType kappa;
};

// Aliases for commonly used types
using DElNinoEq = ElNinoEquation<double>;
using DElNinoSetup = capd::ddeshelper::NonrigorousHelper<DElNinoEq, 1>;
using DElNinoSolution = DElNinoSetup::Solution;
using DElNinoSolver = DElNinoSetup::Solver;
using DElNinoParamType = DElNinoSetup::ParamType;
using DElNinoParamsVector = DElNinoSetup::ParamsVector;
using DElNinoC1Solution = DElNinoSetup::C1Solution;

using IElNinoEq = ElNinoEquation<capd::DInterval>;
using IElNinoSetup = capd::ddeshelper::RigorousHelper<IElNinoEq, 1>;
using IElNinoSolution = IElNinoSetup::Solution;
using IElNinoSolver = IElNinoSetup::Solver;
using IElNinoParamType = IElNinoSetup::ParamType;
using IElNinoParamsVector = IElNinoSetup::ParamsVector;

} // namespace enso
