#include "interpolate.h"
#include <stdexcept>
#include <cmath>

namespace numlab {

    double lagrange(const Vector& xi, const Vector& fi, double x)
    {
        const int n = static_cast<int>(xi.size());
        if (n == 0 || fi.size() != xi.size())
            throw std::runtime_error("Niepoprawne rozmiary wektorow");

        double result = 0.0;

        for (int i = 0; i < n; ++i) {
            double term = fi[i];
            for (int j = 0; j < n; ++j)
                if (j != i)
                    term *= (x - xi[j]) / (xi[i] - xi[j]);
            result += term;
        }
        return result;
    }

    Vector newton_coeff(const Vector& xi, const Vector& fi)
    {
        int n = static_cast<int>(xi.size());
        if (n == 0 || fi.size() != xi.size())
            throw std::runtime_error("Niepoprawne rozmiary wektorow");

        Vector a = fi;           
        for (int j = 1; j < n; ++j)
            for (int i = n - 1; i >= j; --i)
                a[i] = (a[i] - a[i - 1]) / (xi[i] - xi[i - j]);
        return a;               
    }

    double newton_eval(const Vector& a, const Vector& xi, double x)
    {
        int n = static_cast<int>(a.size());
        if (n == 0 || xi.size() != a.size())
            throw std::runtime_error("Niepoprawne rozmiary wektorow");

        double result = a.back();
        for (int i = n - 2; i >= 0; --i)
            result = result * (x - xi[i]) + a[i];
        return result;
    }

    double poly_eval(const Vector& a, double x)
    {
        double result = 0.0;
        for (std::size_t i = 0; i < a.size(); ++i)
            result += a[i] * std::pow(x, static_cast<int>(i));
        return result;
    }

    double poly_horner(const Vector& a, double x)
    {
        if (a.empty()) return 0.0;
        double result = a.back();
        for (int i = static_cast<int>(a.size()) - 2; i >= 0; --i)
            result = result * x + a[i];
        return result;
    }

} 
