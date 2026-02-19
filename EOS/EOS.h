#pragma once

#include<cmath>
#include <vector>
#include "../States/StateW.h"
#include "../Phase/Phase.h"
#include "../Components/Components.h"

inline double ro_P_T(double P, double T, const Phase& ph);
inline double ro_e_P(double e, double P, const Phase& ph);
inline double ro_e_T(double e, double T, const Phase& ph);

inline double P_e_T(double e, double T, const Phase& ph);
inline double P_ro_e(double ro, double e, const Phase& ph);
inline double P_ro_T(double ro, double T, const Phase& ph);

inline double e_P_T(double P, double T, const Phase& ph);
inline double e_ro_P(double ro, double P, const Phase& ph);
inline double e_ro_T(double ro, double T, const Phase& ph);

inline double T_ro_e(double ro, double e, const Phase& ph);
inline double T_ro_P(double ro, double P, const Phase& ph);
inline double T_e_P(double e, double P, const Phase& ph);

// можно оптимизировать вложенные функции

inline double ro_P_T(double P, double T, const Phase& ph)
{
    return (P + ph.P0)/(ph.Cv*T*(ph.gamma -1.0));
}

inline double ro_e_P(double e, double P, const Phase& ph)
{
    return (P + ph.gamma*ph.P0)/((e-ph.q)*(ph.gamma - 1.0));
}

inline double ro_e_T(double e, double T, const Phase& ph)
{
    return (P_e_T(e, T, ph) + ph.gamma*ph.P0)/((e-ph.q)*(ph.gamma - 1.0));
}

inline double P_ro_T(double ro, double T, const Phase& ph)
{
    return T*ph.Cv*(ph.gamma - 1.0) - ph.P0;
}

inline double P_ro_e(double ro, double e, const Phase& ph)
{
    return (e-ph.q)*ro*(ph.gamma - 1.0) - ph.gamma*ph.P0;
}

inline double P_e_T(double e, double T, const Phase& ph)
{
    return (e-ph.q)*ro_e_T(e, T, ph)*(ph.gamma - 1.0) - ph.gamma*ph.P0;
}

inline double e_P_ro(double P, double ro, const Phase& ph)
{
    return (P + ph.gamma*ph.P0)/(ro*(ph.gamma - 1.0)) + ph.q;
};

inline double e_P_T(double P, double T, const Phase& ph)
{
    return (P + ph.gamma*ph.P0)/(ro_P_T(P, T, ph)*(ph.gamma - 1.0)) + ph.q;
};

inline double e_ro_T(double ro, double T, const Phase& ph)
{
     return (P_ro_T(ro, T, ph) + ph.gamma*ph.P0)/(ro*(ph.gamma - 1.0)) + ph.q;
};

inline double T_ro_P(double ro, double P, const Phase& ph)
{
    return (P + ph.P0)/(ph.Cv*ro*(ph.gamma -1.0));
};

inline double T_e_P(double e, double P, const Phase& ph)
{
    return (P + ph.P0)/(ph.Cv*ro_e_P(e, P, ph)*(ph.gamma -1.0));
};

inline double T_ro_e(double ro, double e, const Phase& ph)
{
    return (P_ro_e(ro, e, ph) + ph.P0)/(ph.Cv*ro*(ph.gamma -1.0));
};

inline double S_P_T(double P, double T, const Phase& ph)
{
    return ph.Cv * (ph.gamma*std::log(T) - (ph.gamma - 1.0)*std::log(P + ph.P0)) + ph.qs;
};

inline double E(double e, double u, double v)
{
    return e + (u*u + v*v)*0.5;
};

inline double g_ro_e(double ro, double e, const Phase& ph)
{
    double P = P_ro_e(ro, e, ph);
    double T = T_ro_e(ro, e, ph);
    return e + P/ro - T*S_P_T(P, T, ph);
}

inline double g_P_ro(double P, double ro, const Phase& ph)
{
    double e = e_P_ro(P, ro, ph);
    double T = T_ro_e(ro, e, ph);
    return e + P/ro - T*S_P_T(P, T, ph);
}

inline double f_C2(double P, double ro, const Phase& ph)
{
    return (P/(ro*ro) + ((P+ph.gamma*ph.P0)/(ro*ro*(ph.gamma-1)))) / (1.0/(ro*(ph.gamma-1)));
}

inline double f_C(double P, double ro, const Phase& ph)
{
    return std::sqrt((P/(ro*ro) + ((P+ph.gamma*ph.P0)/(ro*ro*(ph.gamma-1)))) / (1.0/(ro*(ph.gamma-1))));
}

inline double f_PI(const StateW& state)
{
    return state.a1*state.P1 + (1.0-state.a1)*state.P2;
}
//?
inline double f_UI_x(const StateW& state)
{
    
    return (state.a1*state.ro1*state.u1 + (1.0-state.a1)*state.ro2*state.u2)
            /(state.a1*state.ro1 + (1.0-state.a1)*state.ro2);
}

inline double f_UI_y(const StateW& state)
{
    return (state.a1*state.ro1*state.v1 + (1.0-state.a1)*state.ro2*state.v2)
            /(state.a1*state.ro1 + (1.0-state.a1)*state.ro2);
}

inline double T_from_P(double P, double m1, double m2, const Components& comp)
{
    return 1.0 / (m1 * comp.p1.Cv*(comp.p1.gamma-1.0)/(P+comp.p1.P0) 
                + m2 * comp.p2.Cv*(comp.p2.gamma-1.0)/(P+comp.p2.P0));
}