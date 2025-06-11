#include <iostream>
#include <iomanip>
#include <cmath>
#include "integrate.h"

/*****************************************************************************
*  ex_integrate.cpp                                                          *
*                                                                            *
*  Demonstracja modułu integrate:                                            *
*    • całkowanie **dowolnej funkcji f(x)** za pomocą:                       *
*        – midpoint (prostokąty środkowe)                                    *
*        – trapezów                                                          *
*        – Simpsona                                                          *
*        – złożonej kwadratury Gaussa-Legendre’a (3-węzłowej)                *
*    • całkowanie **wielomianu P(x)** opisanego współczynnikami              *
*      – dla wielomianu wywołujemy warianty integral_poly_*                  *
* Funkcje są obiektami z biblioteki Functional                               *
*                                                                            *
*  1.  ∫₃⁵  x cos³x dx   — dokładna wartość  ≈ 6.529 689 124 393 44          *
*  2.  ∫₋₂² P(x) dx ,   P(x)=1·x⁴ − 9·x³ − 5·x² + 0·x − 2                    *
*                                                                            *
*                                                                            *
*****************************************************************************/


using namespace numlab;

/* ---------------- 1. funkcja ciągła ------------------------------------- */
double f(double x) { return x * std::pow(std::cos(x), 3); }

/* ---------------- 2. wielomian 4-go stopnia ----------------------------- *
 *   P(x) = a0 + a1 x + a2 x² + a3 x³ + a4 x⁴                               */
const Vector poly = { -2.0, 0.0, -5.0, -9.0, 1.0 };

int main()
{
    std::cout << std::fixed << std::setprecision(10);

    //Całkowanie funkcji f(x) = x * cos^x
    const double a1 = 3.5, b1 = 6.52968912439344, STEP_N = 100;  // n dla metod klasycznych
    const int n = 100;
    const double nG = 3, mG = 10;                    // Gauss: 3 węzły, 10 paneli, Możliwe są do wyboru także 2, 3 i 4 węzły

    std::cout << "Calka f(x)=x*cos^3(x) na [3.5, 6.52968912439344]\n";
    std::cout << "Midpoint            = "
        << integral_midpoint(f, a1, b1, n) << '\n';
    std::cout << "Trapezy             = "
        << integral_trapezoid(f, a1, b1, n) << '\n';
    std::cout << "Simpson             = "
        << integral_simpson(f, a1, b1, n) << '\n';
    std::cout << "Gauss Legendre (3)  = "
        << integral_gauss_legendre(f, a1, b1, nG, mG) << "\n\n";

    const double a2 = -2.0, b2 = 2.0;  // przedzial dla wielomianu
    std::cout << "Calka wielomianu P(x) na [-2,2]\n";
    std::cout << "Poly Midpoint       = "
        << integral_poly_midpoint(poly, a2, b2, n) << '\n';
    std::cout << "Poly Trapezy        = "
        << integral_poly_trapezoid(poly, a2, b2, n) << '\n';
    std::cout << "Poly Simpson        = "
        << integral_poly_simpson(poly, a2, b2, n) << '\n';
    std::cout << "Poly GaussLeg (3)   = "
        << integral_poly_gauss(poly, a2, b2, nG, mG) << '\n';

    /* UWAGA:
       – dla funkcji ciągłych używamy integral_* (przekazujemy std::function)
       – dla wielomianu w bazie {1,x,x²,…} wygodniej skorzystać z aliasów
         integral_poly_*  (wektor współczynników -> wewnętrznie Horner)      */
}