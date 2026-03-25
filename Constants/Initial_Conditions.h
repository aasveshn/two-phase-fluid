#pragma once
#include "Physics_Constants.h"

const double TL0 = 470.002; //549.998
const double PL0 = 2790000  ;
const double a1L0 = 1e-5;
const double ro1L0 = ro_P_T(PL0, TL0, Phase(gamma_water_vapor,
                                       P0_water_vapor,
                                       Cv_water_vapor,
                                       Cp_water_vapor,
                                       q_water_vapor,
                                       qs_water_vapor));
const double u1L0 = 0;
const double v1L0 = 0;
const double P1L0 = PL0;
const double ro2L0 = ro_P_T(PL0, TL0, Phase(gamma_water_liquid,
                                       P0_water_liquid,
                                       Cv_water_liquid,
                                       Cp_water_liquid,
                                       q_water_liquid,
                                       qs_water_liquid));
const double u2L0 = 0;
const double v2L0 = 0;
const double P2L0 = PL0;

const double TR0 = 300;
const double PR0  = 100000;

const double a1R0 = 1-1e-10;
const double ro1R0 = ro_P_T(PR0, TR0, Phase(gamma_water_vapor,
                                       P0_water_vapor,
                                       Cv_water_vapor,
                                       Cp_water_vapor,
                                       q_water_vapor,
                                       qs_water_vapor));
const double u1R0 = 0;
const double v1R0 = 0;
const double P1R0 = PR0;
const double ro2R0 = ro_P_T(PR0, TR0, Phase(gamma_water_liquid,
                                       P0_water_liquid,
                                       Cv_water_liquid,
                                       Cp_water_liquid,
                                       q_water_liquid,
                                       qs_water_liquid));
const double u2R0 = 0;
const double v2R0 = 0;
const double P2R0 = PR0;


