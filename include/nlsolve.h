#pragma once
#include <vector>
#include <functional>

namespace numlab {

    struct IterData { int iter; double x; };

    constexpr double NL_EPS = 1e-10;   
    constexpr int    NL_MAX_ITER = 100;     

    double root_bisection(const std::function<double(double)>& f,
        double a, double b,
        double eps = NL_EPS,
        int    maxIter = NL_MAX_ITER,
        std::vector<IterData>* trace = nullptr);

    double root_secant(const std::function<double(double)>& f,
        double x0, double x1,
        double eps = NL_EPS,
        int    maxIter = NL_MAX_ITER,
        std::vector<IterData>* trace = nullptr);

    double root_regulafalsi(const std::function<double(double)>& f,
        double a, double b,
        double eps = NL_EPS,
        int    maxIter = NL_MAX_ITER,
        std::vector<IterData>* trace = nullptr);

    double root_newton(const std::function<double(double)>& f,
        const std::function<double(double)>& df,
        double x0,
        double eps = NL_EPS,
        int    maxIter = NL_MAX_ITER,
        std::vector<IterData>* trace = nullptr);

} 

