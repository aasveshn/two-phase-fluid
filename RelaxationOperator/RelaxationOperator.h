#pragma once

#include <iostream>
#include <vector>

#include "../States/StateW.h"
#include "../EOS/EOS.h"
#include "../Mesh/Mesh.h"
#include "../Components/Components.h"


class RelaxationOperator{

private:

    Mesh& mesh;
    const Components& phases;
    double P_min;
    double P_max;
    double P_h;
    std::vector<double> Tsat_curve;

public:

    RelaxationOperator(Mesh& msh, const Components& comp);
    void Relax1();
    void Relax2();
    
private:
    void VelocityRelaxation(Cell& cell);
    void VelocityRelaxation(Cell& cell, bool isX);
    void PressureRelaxation(Cell& cell);
    void PressureTemperatureRelaxation(Cell& cell);
    inline double SolvePressure(double m1, double m2, double E, const Phase& p1, const Phase& p2);
    inline double Fdg(double m1, double m2, double E);
    inline double Fdg(double P, double T);
    inline double Fdg_der(double P, double T, const Phase&p1, const Phase& p2);
    double SolveTemperature(double P, double T0, double tol = 1e-2, int maxIter = 1000);
    inline double GetTsat(double P);
    void GibbsFreeEnergyRelaxation(Cell& cell);

    
};