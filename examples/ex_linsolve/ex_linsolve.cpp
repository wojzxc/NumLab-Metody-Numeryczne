#include <iostream>
#include <iomanip>
#include "linsolve.h"

/******************************************************************************
 *  ex_linsolve.cpp                                                          *
 *  Przykład użycia numlab::gaussian_elimination – rozwiązywanie układów     *
 *  liniowych dowolnego rzędu (nxn) metodą eliminacji Gaussa z częściowym    *
 *  wyborem elementu głównego.                                               *
 *                                                                           *
 *  Jak używać funkcji gaussian_elimination():                               *
 *  -----------------------------------------------------------------------  *
 *  1.   #include "linsolve.h"                                               *
 *  2.   Przygotuj                                                           *
 *         - numlab::Matrix  A  – wektor wierszy (vector<vector<double>>)    *
 *         - numlab::Vector  b  – wektor wyrazów wolnych                     *
 *       Obie struktury muszą mieć zgodne rozmiary (A to nxn, b to n).       *
 *  3.   Wywołaj:  auto x = numlab::gaussian_elimination(A,b);               *
 *       - zwróci wektor rozwiązań; w razie macierzy osobliwej rzuca         *
 *         std::runtime_error.                                               *
 *  4.   Koszt obliczeń O(n^3); algorytm pracuje na kopii A, b (oryginały   *
 *       pozostają nienaruszone).                                            *
 *                                                                           *
 *  Poniżej dwa układy:                                                      *
 *    • mały 3×3 – przykład „podręcznikowy”                                  *
 *    • większy 5×5 – dane z laboratorium                  *
 *                                                                           *
 *  Budowanie w Visual Studio:                                               *
 *    – plik należy do projektu, który linkuje bibliotekę NumLab (Project    *
 *      Reference ▸ NumLab). Wtedy wystarczy F5.                             *
 *    -aby skompilować należy ustawić jako projekt startowy                  *
 ******************************************************************************/


void solve_and_print(const numlab::Matrix& A,
    const numlab::Vector& b,
    const char* title)
{
    auto x = numlab::gaussian_elimination(A, b);

    std::cout << title << "\nRozwiazanie:\n";
    for (std::size_t i = 0; i < x.size(); ++i)
        std::cout << "x[" << i << "] = "
        << std::setw(12) << std::setprecision(6) << x[i] << '\n';
    std::cout << '\n';
}

int main()
{
    /* ---------- 3×3 – szybka demonstracja -------------------------------- */
    numlab::Matrix A3 = {
        {  3,  2,  -1},
        {  2, -2,   4},
        { -1,  0.5,-1}
    };
    numlab::Vector b3 = { 1, -2, 0 };

    solve_and_print(A3, b3, "Uklad 3x3");

    /* ---------- 5×5 – dane z ćwiczeń ------------------------------------- */
    numlab::Matrix A5 = {
        {   0, -10,  16,  15,  -3},
        {  -4,   0,  13,  18,  -7},
        {  -6, -15,   0,  10, -11},
        { -16,  -8,  -2,   0,   5},
        {  11,  17, -12,   5,   0}
    };
    numlab::Vector b5 = { -10, 12, -18, 3, -7 };

    solve_and_print(A5, b5, "Uklad 5x5");
}

