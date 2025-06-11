#pragma once
#include <vector>
#include <functional>

namespace numlab {

    using Vector = std::vector<double>;

    double integral_midpoint(const std::function<double(double)>& f,
        double a, double b, int n);

    double integral_trapezoid(const std::function<double(double)>& f,
        double a, double b, int n);

    double integral_simpson(const std::function<double(double)>& f,
        double a, double b, int n);

    double integral_gauss_legendre(const std::function<double(double)>& f,
        double a, double b,
        int nG = 3, int m = 1);


    double integral_poly_midpoint(const Vector& a, double l, double r, int n);
    double integral_poly_trapezoid(const Vector& a, double l, double r, int n);
    double integral_poly_simpson(const Vector& a, double l, double r, int n);
    double integral_poly_gauss(const Vector& a, double l, double r,
        int nG = 3, int m = 1);

} 

