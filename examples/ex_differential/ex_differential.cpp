/*****************************************************************************
*  ex_differential.cpp                                                       *
*                                                                            *
*  Przykład użycia modułu ode: cztery jawne metody całkowania ODE            *
*                                                                            *
*        dT/dt = -alpha * ( T^4 - Tamb^4 )                                   *
*        T(0)  = 3391 K                                                      *
*        alpha = 1.0e-12 ,  Tamb = 300 K                                     *
*                                                                            *
*  • Euler          – step_euler                                             *
*  • Heun           – step_heun                                              *
*  • Mid-Point      – step_midpoint                                          *
*  • Runge–Kutta 4  – step_rk4                                               *
*                                                                            *
*  Funkcja ode_solve() przyjmuje:                                            *
*      y0, t0, tEnd, h,  lambda f(t,y),  wskaźnik do funkcji kroku           *
*  i zwraca wektor (t, y).                                                   *
*                                                                            *
*  Kompilacja: projekt musi linkować bibliotekę NumLab (Project Reference).  *
*****************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <array>
#include "differential.h"

using namespace numlab;

/* ---------- RHS równania ------------------------------------------------ */
//Definiuje pochodną funkcji T(t) w postaci dT/dt = f(t, T).
double rhs(double /*t*/, double T)
{
    constexpr double alpha = 1.0e-12; 
    constexpr double Tamb = 300.0;
    return -alpha * (std::pow(T, 4) - std::pow(Tamb, 4));
}

/* ---------- pomoc: wypisz kilka punktów trajektorii -------------------- */
void print_traj(const std::vector<StatePoint>& v, int every)
{
    std::cout << " t[s]      T[K]\n";
    for (std::size_t i = 0; i < v.size(); i += every)
        std::cout << std::setw(6) << std::fixed << std::setprecision(0) << v[i].t
        << "   " << std::setprecision(2) << v[i].y << '\n';
    std::cout << "---------------------------------\n";
}

/* ---------- uruchom wybraną metodę ------------------------------------- */
template<typename StepFunc>
void run_one(const char* name, StepFunc step, double h)
{
    auto sol = ode_solve(
        /*y0*/ 3391.0,
        /*t0*/ 0.0,
        /*tEnd*/ 600.0,     // 10 minut
        /*h*/ h,
        rhs,
        step); // wybór kroku

    std::cout << "\n" << name << " ,  h = " << h << " s\n";
    print_traj(sol, 12);    // co 12 kroków (czyli 120 s dla h=10)
}

int main()
{
    const std::array<double, 3> steps{ 1.0, 10.0, 18.0 };

    for (double h : steps)
    {
        run_one("Euler         ", step_euler, h);
        run_one("Heun  (RK2)   ", step_heun, h);
        run_one("Mid-Point (RK2)", step_midpoint, h);
        run_one("Runge-Kutta 4 ", step_rk4, h);
    }
}

