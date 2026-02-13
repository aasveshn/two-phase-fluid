#pragma once

struct StateW {
    
    union {
        struct {
            double a1, ro1, u1, v1, P1, ro2, u2, v2, P2;
        };
        double w[9]; 
    };

    StateW() : a1(0.0), ro1(0.0), u1(0.0), v1(0.0), P1(0.0), ro2(0.0), u2(0.0), v2(0.0), P2(0.0) {}

    StateW(double a1, double r1, double u1, double v1, double P1, double r2, double u2, double v2, double P2)
        : a1(a1), ro1(r1), u1(u1), v1(v1), P1(P1), ro2(r2), u2(u2), v2(v2), P2(P2) {}

   
    inline StateW operator+(const StateW& s) const {
        return {a1 + s.a1, ro1 + s.ro1, u1 + s.u1, v1 + s.v1, P1 + s.P1, 
                ro2 + s.ro2, u2 + s.u2, v2 + s.v2, P2 + s.P2};
    }
    inline StateW operator-(const StateW& s) const {
        return {a1 - s.a1, ro1 - s.ro1, u1 - s.u1, v1-s.v1, P1 - s.P1,
                ro2 - s.ro2, u2 - s.u2, v2-s.v2, P2 - s.P2};
    }
    inline StateW operator*(double c) const {
        return {a1 * c, ro1 * c, u1 * c, v1 * c,  P1 * c,
                ro2 * c, u2 * c, v2 * c, P2 * c};
    }

    inline StateW& operator=(const StateW& other) {
       if (this != &other) { 
            a1 = other.a1; ro1 = other.ro1; u1 = other.u1; v1 = other.v1; P1 = other.P1;
            ro2 = other.ro2; u2 = other.u2; v2 = other.v2; P2 = other.P2;
       }
       return *this;
    }

    
    inline double operator[](int i) const { return w[i]; }
    inline double& operator[](int i) { return w[i]; }
};

inline StateW operator*(double c, const StateW& state) {

    return state * c; 
}