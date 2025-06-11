NumLab  –  biblioteka metod numerycznych


┌────────────────────────────────-------------─┐
│ Moduły:                                      │
│  • linsolve.cpp / .h     – układy liniowe    │
│  • integrate.cpp / .h    – całkowanie num.   │
│  • nlsolve.cpp / .h      – równania nielin.  │
│  • ode.cpp / .h          – ODE 1-go rzędu    │
│  • approx.cpp / .h       – LSQ (MSE) poly    │
│  • interpolate.cpp / .h  – Lagrange / Newton │
└───────────────────────────────-------------──┘


1 BUDOWANIE BIBLIOTEKI (tylko .lib / .a)
-----------------------------------------

CMake (Windows / Linux / macOS)
bash
# konfiguracja (Debug domyślnie)
cmake -S . -B build
# budowa
cmake --build build --target NumLab
• Windows MSVC → build/Debug/NumLab.lib
-----------------------------------
Visual Studio 2022 
Otworzyć projekt sln w Visual Studio -> a następnie CTRL+SHIST+B
po zbudowaniu w katalogu Debug pojawi się plik NumLab.lib



2 DOŁĄCZANIE DO WŁASNEGO PROJEKTU
Visual Studio (sugerowane)

utworzyć folder libs a następnie dołączyć folder include z biblioteki oraz plik NumLab.lib

We właściwościach:
C/C++ -> General ->Additional Include Directories → $(ProjektDir)libs\include
Linker -> General -> Additional Library Directories -> $(ProjektDir)libs 
Linker -> input ->Additional Dependencies → NumLab.lib

Aby użyć właściwą funkcje trzeba w pliku dodać #include "plik nagłówkowy z oczekiwaną funkcją"


3 PRZYKŁADOWE UŻYCIE 


#include <iostream>
#include <cmath>
#include "integrate.h"        // moduł => nagłówek

double f(double x){ return std::cos(x); }

int main(){
    double I = numlab::integral_simpson(f, 0.0, 1.0, 400);
    std::cout << "∫₀¹ cos = " << I << '\\n';     // ≈ 0.841470985
}


4 LISTA FUNKCJI (publiczny interfejs)

-linsolve.h
Funkcja	
Vector gaussian_elimination(Matrix A, Vector b)	rozwiązuje 𝐴 𝑥 = 𝑏
Ax=b z pivotem częściowym	std::runtime_error jeśli macierz osobliwa

-integrate.h
Funkcja	
integral_midpoint(f,a,b,n)	prostokąty (środek)
integral_trapezoid(f,a,b,n)	trapézy
integral_simpson(f,a,b,n)	Simpson (n parzyste auto-poprawka)
integral_gauss_legendre(f,a,b,nG,m)	składany Gauss-Legendre (2/3/4 węzły) ◆
integral_poly_*	analogiczne cztery wersje dla wielomianu podanego współczynnikami Vector

-nlsolve.h
Funkcja	Wymagania	Zwraca
root_bisection(f,a,b)	f(a)·f(b)<0	pierwiastek lub NaN ◆
root_regulafalsi(f,a,b)	j.w.	^
root_secant(f,x0,x1)	dowolny start	pierwiastek
root_newton(f,df,x0)	pochodna df	pierwiastek lub x0 jeśli brak zbieżności

-differential.h
Typ 
StatePoint{t,y}	pojedynczy punkt trajektorii
ode_solve(y0,t0,tEnd,h,f,step)	integrator zwraca vector<StatePoint>
step_euler, step_heun, step_midpoint, step_rk4	pojedyncze kroki (przekazywane do ode_solve)

-approx.h
Funkcja
poly_lsq(f,a,b,m,n)	LSQ – zwraca Vector{a0..am} wielomianu stopnia m
poly_horner(coeff,x)	schemat Hornera dla powyższego

-interpolate.h
Funkcja	
lagrange(xi,fi,x)	wartość wielomianu Lagrange’a
newton_coeff(xi,fi)	różnice dzielone a₀..aₙ
newton_eval(a,xi,x)	wartość wielomianu Newtona
poly_horner(coeff,x)	inline – identyczna jak w approx (pojedyncza definicja)

5 KONWENCJE, WYJĄTKI, JEDNOSTKI
Wszystkie funkcje liczbowe pracują na double (64-bit).

Błędy użytkownika → std::invalid_argument (np. a≥b, nG≠2/3/4, h≤0).

Macierz osobliwa (gaussian) → std::runtime_error.

Funkcje zwracające double mogą zwrócić NaN (brak korzenia / poza domeną).

Żadna funkcja nie modyfikuje argumentów wejściowych (operuje na kopiach).

6 AUTOMATYCZNE TESTY ORAZ PRZYKŁADY
Aby "Odpalić" testy lub examples należy otworzyć plik sln projektu w visual studio, następnie skompilować bibliotekę. 
Następnie kliknąć PPM na wybrany projek np ex_differencial i ustawić jako projekt startowy.

7 LICENCJA
Projekt wyłącznie edukacyjny – brak formalnej licencji.
Możesz kopiować i modyfikować na potrzeby zajęć Metody Numeryczne.