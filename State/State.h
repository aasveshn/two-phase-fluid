#pragma once

struct State {
    
    union {
        struct {
            double a1, ro1, u1, P1, ro2, u2, P2;
        };
        double w[7]; 
    };

    State() : a1(0), ro1(0), u1(0), P1(0), ro2(0), u2(0), P2(0) {}

    State(double a1, double r1, double u1, double P1, double r2, double u2, double P2)
        : a1(a1), ro1(r1), u1(u1), P1(P1), ro2(r2), u2(u2), P2(P2) {}

   
    inline State operator+(const State& s) const {
        return {a1 + s.a1, ro1 + s.ro1, u1 + s.u1, P1 + s.P1, ro2 + s.ro2, u2 + s.u2, P2 + s.P2};
    }
    inline State operator-(const State& s) const {
        return {a1 - s.a1, ro1 - s.ro1, u1 - s.u1, P1 - s.P1, ro2 - s.ro2, u2 - s.u2, P2 - s.P2};
    }
    inline State operator*(double c) const {
        return {a1 * c, ro1 * c, u1 * c, P1 * c, ro2 * c, u2 * c, P2 * c};
    }

    
    inline double operator[](int i) const { return w[i]; }
    inline double& operator[](int i) { return w[i]; }
};
