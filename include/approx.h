#pragma once
#include <vector>
#include <functional>

namespace numlab {

    using Vector = std::vector<double>;

    Vector poly_lsq(const std::function<double(double)>& f,
        double a, double b,
        int    m,
        int    n = 200);

    double poly_horner(const Vector& coef, double x);

} 