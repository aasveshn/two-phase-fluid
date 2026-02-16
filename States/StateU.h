#pragma once

#include "StateW.h"
#include "../EOS/EOS.h"
#include "../Components/Components.h"

struct StateU {
    // Используем union для доступа и по именам, и по индексу
    union {
        struct {
            double ar1, aru1, arv1, arE1, ar2, aru2,arv2, arE2;  
        };
        double v[6]; 
    };

    
    StateU() : ar1(0), aru1(0), arv1(0), arE1(0), ar2(0), aru2(0), arv2(0), arE2(0) {}

    
    StateU(double ar1, double aru1, double arv1, double arE1, 
           double ar2, double aru2, double arv2, double arE2)
        : ar1(ar1), aru1(aru1), arv1(arv1), arE1(arE1), ar2(2), aru2(aru2), arv2(arv2), arE2(arE2) {}

    StateU(const StateW& w, const Components& comp)
    {
        ar1 = w.a1*w.ro1;
        aru1 = w.a1*w.ro1*w.u1;
        arv1 = w.a1*w.ro1*w.v1;
        arE1 = w.a1*w.ro1*E(e_P_ro(w.P1, w.ro1, comp.p1), w.u1, w.v1);
        ar2 = (1-w.a1)*w.ro2;
        aru2 = (1-w.a1)*w.ro2*w.u2;
        arv2 = (1-w.a1)*w.ro2*w.v2;
        arE2 =(1-w.a1)*w.ro2*E(e_P_ro(w.P2, w.ro2, comp.p2), w.u2, w.v2);
    }

    inline void zero() {
        ar1 = aru1 =  arv1 = arE1 = ar2 = aru2 = arv2 = arE2 = 0.0;
    }

    
    inline StateU operator+(const StateU& u) const {
        return StateU(
            ar1 + u.ar1, aru1 + u.aru1, arv1 + u.arv1, arE1 + u.arE1,
            ar2 + u.ar2, aru2 + u.aru2, arv2 + u.arv2, arE2 + u.arE2
        );
    }

    
    inline StateU operator-(const StateU& u) const {
        return StateU(
            ar1 - u.ar1, aru1 - u.aru1, arv1 - u.arv1, arE1 - u.arE1,
            ar2 - u.ar2, aru2 - u.aru2, arv2 - u.arv2, arE2 - u.arE2
        );
    }

    
    inline StateU operator*(double c) const {
        return StateU(
            ar1 * c, aru1 * c, arv1 * c, arE1 * c,
            ar2 * c, aru2 * c, arv2 * c, arE2 * c
        );
    }

    
    inline double operator[](int i) const { return v[i]; }
    inline double& operator[](int i) { return v[i]; }
};

inline StateU get_F_X(const StateU& U, const StateW& W, double a1, double a2) {
    StateU F;
    
    F[0] = U[1];                          
    F[1] = U[1] * W.u1 + a1 * W.P1;       
    F[2] = U[1] * W.v1;                   
    F[3] = (U[3] + a1 * W.P1) * W.u1;     

    F[4] = U[5];                         
    F[5] = U[5] * W.u2 + a2 * W.P2;       
    F[6] = U[5] * W.v2;                  
    F[7] = (U[7] + a2 * W.P2) * W.u2;     

    return F;
}


inline StateU get_F_Y(const StateU& U, const StateW& W, double a1, double a2) {
    StateU G;

    G[0] = U[2];                          
    G[1] = U[2] * W.u1;                   
    G[2] = U[2] * W.v1 + a1 * W.P1;       
    G[3] = (U[3] + a1 * W.P1) * W.v1;     

    G[4] = U[6];                         
    G[5] = U[6] * W.u2;                   
    G[6] = U[6] * W.v2 + a2 * W.P2;       
    G[7] = (U[7] + a2 * W.P2) * W.v2;     

    return G;
}