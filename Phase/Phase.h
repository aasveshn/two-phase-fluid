#pragma once

#include <iostream>

class Phase{

public: 
    double gamma;
    double P0;
    double Cv;
    double Cp;
    double q;
    double qs;

    Phase(double gamma, double P0, double Cv,
          double Cp, double q, double qs);


};


