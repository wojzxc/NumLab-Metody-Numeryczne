#include "linsolve.h"
#include <iostream>
#include <cmath>
#include <iomanip>

namespace {

    void print_step(const numlab::Matrix& M)
    {
        for (const auto& row : M) {
            for (double v : row)
                std::cout << std::setw(10) << std::fixed << std::setprecision(4) << v << ' ';
            std::cout << '\n';
        }
        std::cout << "-------------------------------\n";
    }

} 

namespace numlab {

    Vector gaussian_elimination(Matrix A, Vector b, bool verbose)
    {
        const int n = static_cast<int>(A.size());
        if (n == 0 || static_cast<int>(b.size()) != n)
            throw std::runtime_error("Z³e wymiary uk³adu");

        // 1. tworzymy macierz rozszerzon¹ [A | b]
        for (int i = 0; i < n; ++i)
            A[i].push_back(b[i]);

        if (verbose) {
            std::cout << "Macierz rozszerzona [A|b] - przed eliminacj¹:\n";
            print_step(A);
        }

        for (int i = 0; i < n; ++i)
        {
            int maxRow = i;
            for (int k = i + 1; k < n; ++k)
                if (std::fabs(A[k][i]) > std::fabs(A[maxRow][i]))
                    maxRow = k;
            std::swap(A[i], A[maxRow]);

            if (std::fabs(A[i][i]) < 1e-12)
                throw std::runtime_error("Macierz osobliwa — brak rozwi¹zania");

            for (int k = i + 1; k < n; ++k)
            {
                double factor = A[k][i] / A[i][i];
                for (int j = i; j <= n; ++j)
                    A[k][j] -= factor * A[i][j];

                if (verbose) {
                    std::cout << "Po wyzerowaniu elementu w wierszu "
                        << k << ", kolumnie " << i << ":\n";
                    print_step(A);
                }
            }
        }

        Vector x(n);
        for (int i = n - 1; i >= 0; --i)
        {
            x[i] = A[i][n];
            for (int j = i + 1; j < n; ++j)
                x[i] -= A[i][j] * x[j];
            x[i] /= A[i][i];
        }

        return x;
    }

} 
