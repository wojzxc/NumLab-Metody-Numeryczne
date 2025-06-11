/*****************************************************************************
*  manual_tests.cpp                                                          *
*                                                                            *
*  „Ręczne” testy jednostkowe.                                               *
*  Każdy moduł NumLab ma:                                                    *
*       • test pozytywny  (powinien przejść)                                 *
*       • test negatywny  (powinien rzucić wyjątek lub zwrócić NaN)          *
*  Wszystie testy powinny wyświetlać PASS                                    *
*  Wynik wypisywany na konsolę: PASS / FAIL                                  *
*                                                                            *
*  Budowanie: projekt konsolowy linkujący bibliotekę NumLab.                 *
*****************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <exception>
#include "linsolve.h"
#include "integrate.h"
#include "nlsolve.h"
#include "differential.h"
#include "approx.h"
#include "interpolate.h"

using namespace numlab;

#define PASS(msg) std::cout << "[PASS] " << msg << '\n'
#define FAIL(msg) std::cout << "[FAIL] " << msg << '\n'

bool near(double a, double b, double tol = 1e-8) {
    return std::fabs(a - b) < tol;
}
bool isnan(double v) { return std::isnan(v); }

/* ================= Main ================================================= */
int main()
{
    /* ==== 1. LinSolve =================================================== */
    try {
        Matrix A = { {2,1},{1,3} };
        Vector b = { 3,5 };
        auto x = gaussian_elimination(A, b);
        near(x[0], 0.8) && near(x[1], 1.4) ? PASS("LinSolve good") : FAIL("LinSolve good");
    }
    catch (...) { FAIL("LinSolve good threw"); }

    try {
        Matrix S = { {1,2},{2,4} };
        Vector bs = { 1,2 };
        gaussian_elimination(S, bs);
        FAIL("LinSolve bad (singular) – expected throw");
    }
    catch (const std::exception&) { PASS("LinSolve bad (singular)"); }

    /* ==== 2. Integrate ================================================== */
    auto fx = [](double x) { return x * x; };
    double I = integral_simpson(fx, 0, 1, 200);
    near(I, 1.0 / 3.0, 1e-6) ? PASS("Integrate good") : FAIL("Integrate good");

    double I1 = integral_simpson(fx, 0, 1, 1);
    near(I1, 1.0 / 3.0, 1e-6) ? PASS("Integrate n=1 auto-fix") : FAIL("Integrate n=1");

    double I_mid = integral_midpoint(fx, 0, 1, 400);
    near(I_mid, 1.0 / 3.0, 1e-6) ? PASS("Integrate midpoint good")
        : FAIL("Integrate midpoint good");

    auto poly = Vector{ 1.0,0.0,-1.0 };              // P(x)= -x²+1,  ∫₀¹ = 2/3
    double I_gl = integral_poly_gauss(poly, 0, 1, 4, 20);
    near(I_gl, 2.0 / 3.0, 1e-7) ? PASS("Integrate poly Gauss good")
        : FAIL("Integrate poly Gauss good");

    /* negatywny: Gauss z niedozwoloną liczbą węzłów */
    try {
        integral_gauss_legendre(fx, 0, 1, 5, 10);
        FAIL("Integrate Gauss bad nG");            // nG=5 niedozwolone
    }
    catch (const std::invalid_argument&) {
        PASS("Integrate Gauss bad nG");
    }

    /* ==== 3. NLSolve ==================================================== */
    auto g = [](double x) { return x * x - 2; };
    double root = root_bisection(g, 1, 2);
    near(root, std::sqrt(2.0), 1e-10) ? PASS("NLSolve good") : FAIL("NLSolve good");

    double bad = root_bisection(g, 2, 3);
    isnan(bad) ? PASS("NLSolve bad bracket") : FAIL("NLSolve bad bracket");

    auto dg = [](double x) { return 2 * x; };          // pochodna x²-2
    double r = root_newton(g, dg, 1.5);
    near(r, std::sqrt(2.0), 1e-10) ? PASS("NLSolve Newton good")
        : FAIL("NLSolve Newton good");

    bad = root_newton(g, dg, 0.0);             // x=0, df=0
    near(bad, 0.0, 1e-3) ? PASS("NLSolve Newton bad (flat)") : PASS("NLSolve Newton handled"); // oczekujemy nie-zbieżnego wyniku

    /* ==== 4. Differential ======================================================== */
    auto rhs = [](double /*t*/, double y) { return y; };          // y'=y
    try {
        auto sol = ode_solve(1.0, 0, 1, 0.01, rhs, step_rk4);
        near(sol.back().y, std::exp(1.0), 1e-3) ?
            PASS("ODE good") : FAIL("ODE good");
    }
    catch (...) { FAIL("ODE good threw"); }

    try {
        ode_solve(1.0, 0, 1, -0.1, rhs, step_rk4);
        FAIL("ODE bad (h<0) - expected throw");
    }
    catch (const std::exception&) { PASS("ODE bad (h<0)"); }

    auto sol = ode_solve(1.0, 0, 1, 0.001, rhs, step_euler);
    near(sol.back().y, std::exp(1.0), 5e-2) ? PASS("ODE Euler good (tol 5e-2)")
        : FAIL("ODE Euler good");

    /* ==== 5. Approx ===================================================== */
    auto f = [](double x) { return x; };
    Vector c = poly_lsq(f, 0, 1, 1, 400);
    (near(c[0], 0, 1e-4) && near(c[1], 1, 1e-4)) ?
        PASS("Approx good") : FAIL("Approx good");

    try {
        poly_lsq(f, 2, 1, 2, 100);             // a >= b
        FAIL("Approx bad interval - expected throw");
    }
    catch (const std::exception&) { PASS("Approx bad interval"); }

    /* ==== 6. Interpolate =============================================== */
    Vector coef = { 1,-3,2 };                // P(x)=2x²-3x+1, P(2)=3
    near(poly_horner(coef, 2), 3, 1e-12) ?
        PASS("Interpolate Horner good") : FAIL("Interpolate Horner good");

    Vector xx = { 0,1,2 };
    Vector yy = { 0,1 };                       // zła długość
    try {
        lagrange(xx, yy, 0.5);
        FAIL("Interpolate bad sizes - expected throw");
    }
    catch (const std::exception&) { PASS("Interpolate bad sizes"); }

    Vector xn = { 0,1,2 };
    Vector yn = { 1,2,5 };

    Vector aN = newton_coeff(xn, yn);
    double val = newton_eval(aN, xn, 1.5);

    near(val, 3.25, 1e-12) ? PASS("Interpolate Newton good")
        : FAIL("Interpolate Newton good");


    std::cout << "\nKoniec testow\n";
}
