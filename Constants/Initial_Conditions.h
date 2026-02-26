#pragma once
#include "Physics_Constants.h"

const double TL0 = 503.0;
const double PL0 = 2790000;
const double a1L0 = 1e-8;
const double ro1L0 = ro_P_T(PL0, TL0, Phase(gamma_dodecane_vapor,
                                       P0_dodecane_vapor,
                                       Cv_dodecane_vapor,
                                       Cp_dodecane_vapor,
                                       q_dodecane_vapor,
                                       qs_dodecane_vapor));
const double u1L0 = 0;
const double v1L0 = 0;
const double P1L0 = 2790000;
const double ro2L0 = ro_P_T(PL0, TL0, Phase(gamma_dodecane_liquid,
                                       P0_dodecane_liquid,
                                       Cv_dodecane_liquid,
                                       Cp_dodecane_liquid,
                                       q_dodecane_liquid,
                                       qs_dodecane_liquid));
const double u2L0 = 0;
const double v2L0 = 0;
const double P2L0 = 2790000;

const double TR0 = 300;
const double PR0  = 100000;

const double a1R0 = 1-1e-8;
const double ro1R0 = ro_P_T(PR0, TR0, Phase(gamma_dodecane_vapor,
                                       P0_dodecane_vapor,
                                       Cv_dodecane_vapor,
                                       Cp_dodecane_vapor,
                                       q_dodecane_vapor,
                                       qs_dodecane_vapor));
const double u1R0 = 0;
const double v1R0 = 0;
const double P1R0 = PR0;
const double ro2R0 = ro_P_T(PR0, TR0, Phase(gamma_dodecane_liquid,
                                       P0_dodecane_liquid,
                                       Cv_dodecane_liquid,
                                       Cp_dodecane_liquid,
                                       q_dodecane_liquid,
                                       qs_dodecane_liquid));
const double u2R0 = 0;
const double v2R0 = 0;
const double P2R0 = PR0;