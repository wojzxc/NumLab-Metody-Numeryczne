#include "approx.h"
#include "linsolve.h"     
#include "integrate.h"
#include <cmath>
#include <stdexcept>

namespace numlab {

    static double simpson(const std::function<double(double)>& g,
        double a, double b, int n)
    {
        if (n % 2) ++n;                       
        double h = (b - a) / n;
        double sum = g(a) + g(b);
        for (int k = 1; k < n; ++k) {
            double x = a + k * h;
            sum += (k % 2 ? 4.0 : 2.0) * g(x);
        }
        return sum * h / 3.0;
    }

    Vector poly_lsq(const std::function<double(double)>& f,
        double a, double b, int m, int n)
    {
        if (m < 0)            throw std::invalid_argument("stopień < 0");
        if (a >= b)           throw std::invalid_argument("a ≥ b");
        if (n < 2)            throw std::invalid_argument("n < 2");

        Matrix A(m + 1, Vector(m + 1));
        Vector B(m + 1);

        for (int i = 0; i <= m; ++i)
            for (int j = 0; j <= m; ++j)
                A[i][j] = simpson([=](double x) { return std::pow(x, i + j); },
                    a, b, n);

        for (int i = 0; i <= m; ++i)
            B[i] = simpson([&](double x) { return f(x) * std::pow(x, i); },
                a, b, n);

        return gaussian_elimination(A, B);
    }


} 
