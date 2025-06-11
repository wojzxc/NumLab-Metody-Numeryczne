#pragma once
#include <vector>
#include <functional>

namespace numlab {

    struct StatePoint { double t, y; };

    using OdeRHS = std::function<double(double, double)>;

    using OdeStep = std::function<double(double, double, double, const OdeRHS&)>;

    double step_euler(double y, double t, double h, const OdeRHS& f);
    double step_heun(double y, double t, double h, const OdeRHS& f);  
    double step_midpoint(double y, double t, double h, const OdeRHS& f);  
    double step_rk4(double y, double t, double h, const OdeRHS& f);  

    std::vector<StatePoint>
        ode_solve(double y0,            
            double t0, double tEnd,
            double h,
            const OdeRHS& f,
            const OdeStep& step = step_rk4);

} 

