#include "nlsolve.h"
#include <cmath>
#include <limits>

namespace numlab {

    double root_bisection(const std::function<double(double)>& f,
        double a, double b, double eps, int maxIter,
        std::vector<IterData>* trace)
    {
        double fa = f(a), fb = f(b);
        if (!std::isfinite(fa) || !std::isfinite(fb) || fa * fb > 0) return std::numeric_limits<double>::quiet_NaN();

        for (int k = 1; k <= maxIter; ++k) {
            double c = 0.5 * (a + b), fc = f(c);
            if (trace) trace->push_back({ k,c });
            if (std::fabs(fc) < eps || 0.5 * (b - a) < eps) return c;
            (fa * fc < 0) ? (b = c, fb = fc) : (a = c, fa = fc);
        }
        return 0.5 * (a + b);
    }

    double root_secant(const std::function<double(double)>& f,
        double x0, double x1, double eps, int maxIter,
        std::vector<IterData>* trace)
    {
        double f0 = f(x0), f1 = f(x1);
        for (int k = 1; k <= maxIter; ++k) {
            if (std::fabs(f1 - f0) < 1e-14) break;
            double x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
            if (trace) trace->push_back({ k,x2 });
            if (std::fabs(x2 - x1) < eps) return x2;
            x0 = x1; f0 = f1; x1 = x2; f1 = f(x1);
        }
        return x1;
    }

    double root_regulafalsi(const std::function<double(double)>& f,
        double a, double b, double eps, int maxIter,
        std::vector<IterData>* trace)
    {
        double fa = f(a), fb = f(b); if (!std::isfinite(fa) || !std::isfinite(fb) || fa * fb > 0) return NAN;
        double x = a;
        for (int k = 1; k <= maxIter; ++k) {
            x = (a * fb - b * fa) / (fb - fa);
            double fx = f(x);
            if (trace) trace->push_back({ k,x });
            if (std::fabs(fx) < eps) return x;
            (fa * fx < 0) ? (b = x, fb = fx) : (a = x, fa = fx);
        }
        return x;
    }

    double root_newton(const std::function<double(double)>& f,
        const std::function<double(double)>& df,
        double x0, double eps, int maxIter,
        std::vector<IterData>* trace)
    {
        for (int k = 1; k <= maxIter; ++k) {
            double dfx = df(x0); if (std::fabs(dfx) < 1e-14) break;
            double x1 = x0 - f(x0) / dfx;
            if (trace) trace->push_back({ k,x1 });
            if (std::fabs(x1 - x0) < eps) return x1;
            x0 = x1;
        }
        return x0;
    }

} 
