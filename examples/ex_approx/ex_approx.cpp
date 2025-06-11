#include <iostream>
#include <iomanip>
#include <cmath>
#include "approx.h"
#include "interpolate.h"

/*****************************************************************************
*  ex_approx.cpp                                                             *
*                                                                            *
*  Aproksymacja średniokwadratowa funkcji                                    *
*        f(x) = sin(x)    na przedziale  [0 , pi/2]                          *
*  wielomianem stopnia m (tutaj m = 3).                                      *
*                                                                            *
*  KROKI:                                                                    *
*   1.  definiujemy lambdę f(x);                                             *
*   2.  wywołujemy                                                           *
*         Vector c = numlab::poly_lsq(f, a, b, m, n);                        *
*       gdzie n = liczba podprzedziałów w regule Simpsona;                   *
*   3.  korzystamy z numlab::poly_horner(c,x) żeby szybko policzyć P(x).     *
*                                                                            *
*  Plik pokazuje:                                                            *
*   • jak uzyskać współczynniki a0..am,                                      *
*   • jak obliczać przybliżenie i błąd rms na siatce kontrolnej.             *
*****************************************************************************/ 


using namespace numlab;

#define M_PI_2 1.57079632679489661923 // pi/2

int main()
{
    // ---- definicja funkcji i parametrów aproksymacji ---- 
    //Obiekt z biblioteki functional, który przyjmuje x i zwraca sin(x).
    auto f = [](double x) { return std::sin(x); };

    const double a = 0.0;
    const double b = M_PI_2;        //  pi/2
    const int    m = 2;             // stopień wielomianu
    const int    n = 800;           // podprzedziały Simpsona

    // ---- współczynniki wielomianu ---- 
    Vector coef = poly_lsq(f, a, b, m, n);

    std::cout << "Wielomian stopnia " << m << " dla sin(x) na [0,pi/2]\n";
    for (int i = 0; i <= m; ++i)
        std::cout << "a" << i << " = " << std::setprecision(10) << coef[i] << '\n';

    // ---- test w kilku punktach ---- 
    std::cout << "\nPorownanie f(x) i P(x):\n"
        << "  x        f(x)         P(x)        abs err\n";

    std::cout << std::fixed << std::setprecision(6);
    for (double x : {0.0, 0.3, 0.6, 1.0})
    {
        double fx = f(x);
        double px = poly_horner(coef, x);
        std::cout << std::setw(5) << x
            << "   " << std::setw(10) << fx
            << "   " << std::setw(10) << px
            << "   " << std::fabs(fx - px) << '\n';
    }

    // ---- RMS-blad na siatce 100 punktow kontrolnych ---- 
    double rms = 0.0;
    for (int k = 0; k <= 100; ++k)
    {
        double x = a + (b - a) * k / 100.0;
        double e = f(x) - poly_horner(coef, x);
        rms += e * e;
    }
    rms = std::sqrt(rms / 101.0);
    std::cout << "\nRMS blad (101 punktow kontrolnych) = " << rms << '\n';

    /* PS:
       poly_lsq  – zwraca vector<double> z a0..am
       poly_horner – szybkie P(x) dla dowolnego x
       Zmieniając m lub n, można poprawiać / badać dokładność */
}

