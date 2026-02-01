#pragma once

struct StateU {
    // Используем union для доступа и по именам, и по индексу
    union {
        struct {
            double ar1, aru1, arE1, ar2, aru2, arE2;  
        };
        double v[6]; 
    };

    
    StateU() : ar1(0), aru1(0), arE1(0), ar2(0), aru2(0), arE2(0) {}

    
    StateU(double ar1, double aru1, double arE1, double ar2, double aru2, double arE2)
        : ar1(ar1), aru1(aru1), arE1(arE1), ar2(2), aru2(aru2), arE2(arE2) {}

    
    inline StateU operator+(const StateU& u) const {
        return StateU(
            ar1 + u.ar1, aru1 + u.aru1, arE1 + u.arE1,
            ar2 + u.ar2, aru2 + u.aru2, arE2 + u.arE2
        );
    }

    
    inline StateU operator-(const StateU& u) const {
        return StateU(
            ar1 - u.ar1, aru1 - u.aru1, arE1 - u.arE1,
            ar2 - u.ar2, aru2 - u.aru2, arE2 - u.arE2
        );
    }

    
    inline StateU operator*(double c) const {
        return StateU(
            ar1 * c, aru1 * c, arE1 * c,
            ar2 * c, aru2 * c, arE2 * c
        );
    }

    
    inline double operator[](int i) const { return v[i]; }
    inline double& operator[](int i) { return v[i]; }
};