#pragma once

#include<cmath>
#include <vector>
#include "../States/StateW.h"
#include "../Phase/Phase.h"
#include "../Components/Components.h"

//добавить "прегрузки"
inline double e_P_ro(double P, double ro, const Phase& ph)
{
    return (P + ph.gamma*ph.P0)/(ro*(ph.gamma - 1.0)) + ph.q;
};
//добавить "прегрузки"
inline double T_P_ro(double P, double ro, const Phase& ph)
{
    return (P + ph.P0)/(ph.Cv*ro*(ph.gamma -1.0));
};
//добавить "прегрузки"
inline double S_P_T(double P, double T, const Phase& ph)
{
    return ph.Cv * (ph.gamma*std::log(T) - (ph.gamma - 1.0)*std::log(P + ph.P0)) + ph.qs;
};
//добавить "прегрузки"
inline double E(double e, double u)
{
    return e + (u*u)*0.5;
};

