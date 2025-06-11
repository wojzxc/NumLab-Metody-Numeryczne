#include "differential.h"

namespace numlab {

    double step_euler(double y, double t, double h, const OdeRHS& f)
    {
        return y + h * f(t, y);
    }
    double step_heun(double y, double t, double h, const OdeRHS& f)
    {
        double k1 = f(t, y);
        double k2 = f(t + h, y + h * k1);
        return y + 0.5 * h * (k1 + k2);
    }
    double step_midpoint(double y, double t, double h, const OdeRHS& f)
    {
        double k1 = f(t, y);
        double k2 = f(t + 0.5 * h, y + 0.5 * h * k1);
        return y + h * k2;
    }
    double step_rk4(double y, double t, double h, const OdeRHS& f)
    {
        double k1 = f(t, y);
        double k2 = f(t + 0.5 * h, y + 0.5 * h * k1);
        double k3 = f(t + 0.5 * h, y + 0.5 * h * k2);
        double k4 = f(t + h, y + h * k3);
        return y + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
    }

    std::vector<StatePoint>
        ode_solve(double y0, double t0, double tEnd, double h,
            const OdeRHS& f, const OdeStep& step)
    {
        int N = static_cast<int>((tEnd - t0) / h + 0.5);
        std::vector<StatePoint> traj;
        traj.reserve(N + 1);

        double t = t0, y = y0;
        for (int i = 0; i <= N; ++i) {
            traj.push_back({ t,y });
            y = step(y, t, h, f);
            t += h;
        }
        if (t < tEnd - 1e-12)                 
            traj.push_back({ tEnd, step(y,t, tEnd - t, f) });
        return traj;
    }

} 
