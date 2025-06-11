#include "integrate.h"
#include <cmath>
#include <stdexcept>

namespace numlab {

    static double horner(const Vector& c, double x)
    {
        if (c.empty()) return 0.0;
        double res = c.back();
        for (int i = static_cast<int>(c.size()) - 2; i >= 0; --i)
            res = res * x + c[i];
        return res;
    }

    double integral_midpoint(const std::function<double(double)>& f,
        double a, double b, int n)
    {
        if (n <= 0) throw std::invalid_argument("n <= 0");
        double h = (b - a) / n, sum = 0.0;
        for (int i = 0; i < n; ++i)
            sum += f(a + (i + 0.5) * h);
        return sum * h;
    }

    double integral_trapezoid(const std::function<double(double)>& f,
        double a, double b, int n)
    {
        if (n <= 0) throw std::invalid_argument("n <= 0");
        double h = (b - a) / n;
        double sum = 0.5 * (f(a) + f(b));
        for (int i = 1; i < n; ++i)
            sum += f(a + i * h);
        return sum * h;
    }

    double integral_simpson(const std::function<double(double)>& f,
        double a, double b, int n)
    {
        if (n % 2) ++n;                       
        double h = (b - a) / n, sum = f(a) + f(b);
        for (int i = 1; i < n; ++i) {
            double x = a + i * h;
            sum += (i % 2 ? 4.0 : 2.0) * f(x);
        }
        return sum * h / 3.0;
    }

    static const double X2[] = { -1 / std::sqrt(3.0),  1 / std::sqrt(3.0) };
    static const double W2[] = { 1.0, 1.0 };
    static const double X3[] = { -std::sqrt(3.0 / 5.0), 0.0,  std::sqrt(3.0 / 5.0) };
    static const double W3[] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
    static const double X4[] = { -0.861136, -0.339981, 0.339981, 0.861136 };
    static const double W4[] = { 0.347855, 0.652145, 0.652145, 0.347855 };

    double integral_gauss_legendre(const std::function<double(double)>& f,
        double a, double b, int nG, int m)
    {
        const double* X, * W; int k;
        if (nG == 2) { X = X2; W = W2; k = 2; }
        else if (nG == 3) { X = X3; W = W3; k = 3; }
        else if (nG == 4) { X = X4; W = W4; k = 4; }
        else throw std::invalid_argument("nG musi być 2,3 lub 4");

        if (m <= 0) throw std::invalid_argument("m <= 0");
        double h = (b - a) / m, sum = 0.0;

        for (int j = 0; j < m; ++j) {
            double aj = a + j * h, bj = aj + h;
            double mid = 0.5 * (aj + bj), half = 0.5 * (bj - aj);
            for (int i = 0; i < k; ++i)
                sum += W[i] * f(mid + half * X[i]);
        }
        return sum * (h / 2.0);
    }

    static std::function<double(double)> make_poly(const Vector& a)
    {
        return [&a](double x) { return horner(a, x); };
    }

    double integral_poly_midpoint(const Vector& a, double l, double r, int n) {
        return integral_midpoint(make_poly(a), l, r, n);
    }
    double integral_poly_trapezoid(const Vector& a, double l, double r, int n) {
        return integral_trapezoid(make_poly(a), l, r, n);
    }
    double integral_poly_simpson(const Vector& a, double l, double r, int n) {
        return integral_simpson(make_poly(a), l, r, n);
    }
    double integral_poly_gauss(const Vector& a, double l, double r,
        int nG, int m) {
        return integral_gauss_legendre(make_poly(a), l, r, nG, m);
    }

} 
