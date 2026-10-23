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
inline T get_pi_limit() {
    return T::pi();
}

template<>
inline double get_pi_limit() {
    return 3.14159265358979323846;
}

template<typename T>
inline T signum(T const& x) {
    if (x > 0) {
        return T(1.0);
    } else if (x < 0) {
        return T(-1);
    } else {
        return T(0);
    }
}

template<>
inline capd::interval signum(capd::interval const& x) {
    if (x > 0) {
        return capd::interval {1.0};
    } else if (x < 0) {
        return capd::interval {-1.0};
    } else if (x.rightBound() != 0 && x.leftBound() != 0) {
        return capd::interval {-1.0, 1.0};
    } else if (x.rightBound() == 0 && x.leftBound() != 0) {
        return capd::interval {-1.0, 0.0};
    } else if (x.rightBound() != 0 && x.leftBound() == 0) {
        return capd::interval {0.0, 1.0};
    } else {
        return capd::interval {0.0};
    }
}

template<class Scalar, class Param = Scalar>
class ElNinoLimitEquation {
public:
    using ScalarType = Scalar;
    using size_type = unsigned int;
    typedef ScalarType RealType;
    using ParamType = Param;
    using VectorType = capd::vectalg::Vector<ScalarType, 0>;
    using ParamsVectorType = capd::vectalg::Vector<ParamType, 0>;

    ElNinoLimitEquation() : alpha(1), beta(1) {}

    ElNinoLimitEquation(ParamType alpha, ParamType beta, ParamType kappa) :
        alpha(alpha),
        beta(beta) {}

    ElNinoLimitEquation(capd::vectalg::Vector<ParamType, 0> const& params) :
        alpha(params[0]),
        beta(params[1]) {}

    static size_type imageDimension() {
        return 1;
    };

    static size_type dimension() {
        return 2;
    }

    static size_type getParamsCount() {
        return 2;
    }

    template<typename Real, typename InVector, typename OutVector>
    void operator()(Real const& t, InVector const x, OutVector& fx) const {
        // a hack to make depreciation warning about implicit copy construction
        // go away
        // instead of 'auto Ttau = x[1]' we create an empty object and then copy
        // with operator=
        typename InVector::ScalarType Ttau;
        Ttau = x[1];
        auto const pi = get_pi_limit<ScalarType>();
        fx[0] = -alpha * signum(Ttau) + beta * cos(2 * pi * t);
    }

    static std::string show() {
        return "El Nino Limit equation as autonomous system: $T'(t) = \\beta \\cos(2 \\pi s) - \\alpha \\sgn( T(t-\\tau))$.";
    }

protected:
    ParamType alpha;
    ParamType beta;
    ParamType kappa;
};

// Aliases for commonly used types
using DElNinoLimitEq = ElNinoLimitEquation<double>;
using DElNinoLimitSetup =
    capd::ddeshelper::NonrigorousHelper<DElNinoLimitEq, 1>;
using DElNinoLimitSolution = DElNinoLimitSetup::Solution;
using DElNinoLimitSolver = DElNinoLimitSetup::Solver;
using DElNinoLimitParamType = DElNinoLimitSetup::ParamType;
using DElNinoLimitParamsVector = DElNinoLimitSetup::ParamsVector;
using DElNinoLimitC1Solution = DElNinoLimitSetup::C1Solution;

using IElNinoLimitEq = ElNinoLimitEquation<capd::DInterval>;
using IElNinoLimitSetup = capd::ddeshelper::RigorousHelper<IElNinoLimitEq, 1>;
using IElNinoLimitSolution = IElNinoLimitSetup::Solution;
using IElNinoLimitSolver = IElNinoLimitSetup::Solver;
using IElNinoLimitParamType = IElNinoLimitSetup::ParamType;
using IElNinoLimitParamsVector = IElNinoLimitSetup::ParamsVector;

} // namespace enso
