#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "nlsolve.h"

/*****************************************************************************
*  ex_nlsolve.cpp                                                            *
*                                                                            *
*  Pokazuje użycie modułu nlsolve do wyszukania **wszystkich** pierwiastków  *
*  w przedziale [-3 , 3] dla dwóch niezależnych równań:                      *
*                                                                            *
*    (1)  log(x + 1) - x^3          = 0      (domena: x > -1)                *
*    (2)  cosh(x) - |x| - 0.55       = 0      (domena: R)                    *
*                                                                            *
*  Metoda:                                                                   *
*    • dzielimy przedział co 0.1,                                            *
*    • sprawdzamy zmianę znaku f(a)*f(b) < 0,                                *
*    • jeżeli wystąpi – uruchamiamy bisekcję z modułu nlsolve,               *
*    • korzenie zbieramy do wektora unikając duplikatów (tol = 1e-9).        *  
*                                                                            *              
*   Funkcje są obiektami z biblioteki functional                             *
*  Uwaga na domenę równania (1): przedział [-3, -1] pomijamy, bo log(x+1)    *
*  jest tam niezdefiniowany.                                                 *
*****************************************************************************/


using namespace numlab;

//W przypadku metody Newtona potrzebujemy funkcji f(x) i jej pochodnej df(x).
/* ---------------- równanie 1 ------------------------------------------- */
double f1(double x) { return std::log(x + 1.0) - std::pow(x, 3); }
double df1(double x) { return 1.0 / (x + 1.0) - 3.0 * std::pow(x, 2); }

/* ---------------- równanie 2 ------------------------------------------- */
double f2(double x) { return std::cosh(x) - std::fabs(x) - 0.55; }
double df2(double x) { return std::sinh(x) - (x >= 0 ? 1 : -1); }

/* ----------------------------------------------------------------------- */
constexpr double LEFT = -3.0;
constexpr double RIGHT = 3.0;
constexpr double STEP = 0.1; //Wybieramy możliwie mały krok, aby znaleźć wszystkie pierwiastki

/* unikaj dubletów */
static void add_unique(std::vector<double>& v, double r, double tol = 1e-9) {
    for (double q : v) if (std::fabs(q - r) < tol) return;
    v.push_back(r);
}

/* uniwersalna procedura: skan + 4 funkcje */
template<typename F, typename DF>
void run_all_methods(const char* title, F f, DF df)
{
    std::cout << "\n=== " << title << " ===\n"
        << std::setw(14) << "bisekcja"
        << std::setw(14) << "falsi"
        << std::setw(14) << "sieczne "
        << "    Newton\n";

    double a = LEFT, b = LEFT + STEP;
    while (b <= RIGHT + 1e-12) {
        double fa = f(a), fb = f(b);

        /* pomijamy podprzedziały, gdzie f jest NaN (np. x<=-1 dla f1) */
        if (std::isfinite(fa) && std::isfinite(fb) && fa * fb < 0.0)
        {
            double r_bis = root_bisection(f, a, b);
            double r_falsi = root_regulafalsi(f, a, b);
            double r_sec = root_secant(f, a, b);
            double r_new = root_newton(f, df, 0.5 * (a + b));

            std::cout <<  std::setw(14) << r_bis << ' '
                << std::setw(14) << r_falsi << ' '
                << std::setw(14) << r_sec << ' '
                << r_new << '\n';
        }
        a = b; b += STEP;
    }
}

int main()
{
    std::cout << std::fixed << std::setprecision(10);

    run_all_methods("log(x+1) - x^3", f1, df1);
    run_all_methods("cosh(x) - |x| - 0.55", f2, df2);
}
