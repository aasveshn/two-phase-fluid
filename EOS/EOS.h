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

inline double E(double e, double u)
{
    return e + (u*u)*0.5;
};

inline double g(double ro, double e, const Phase& ph)
{
    double P = P_ro_e(ro, e, ph);
    double T = T_ro_e(ro, e, ph);
    return e + P/ro - T*S_P_T(P, T, ph);
}