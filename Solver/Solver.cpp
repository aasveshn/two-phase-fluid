#include "Solver.h"

Solver::Solver(Mesh& msh):mesh(msh),
                          phases(Phase(gamma_water_vapor,
                                       P0_water_vapor,
                                       Cv_water_vapor,
                                       Cp_water_vapor,
                                       q_water_vapor,
                                       qs_water_vapor),
                                 Phase(gamma_water_liquid,
                                       P0_water_liquid,
                                       Cv_water_liquid,
                                       Cp_water_liquid,
                                       q_water_liquid,
                                       qs_water_liquid)), 
                           HyperbolicOp(msh, phases) {}

double  Solver::compute_dt()
{
    std::vector<double> EightValues(7);
    double h = std::max(dx, dy);
    double MaxEihgtValue = -1;
    for(int i = 0; i < mesh.Cells.size(); ++i)
    {
        StateW& W = mesh.Cells[i].W;
        double C1 = f_C(W.P1, W.ro1, phases.p1);
        double C2 = f_C(W.P2, W.ro2, phases.p2);
        double UI = f_UI(W);

        EightValues[0] = UI;
        EightValues[1] = W.u1 - C1;
        EightValues[2] = W.u1;
        EightValues[3] = W.u1 + C1;
        EightValues[4] = W.u2 - C2;
        EightValues[5] = W.u2;
        EightValues[6] = W.u2 + C2;

        for(int j = 0; j < 7; ++j)
            MaxEihgtValue = std::max(MaxEihgtValue, std::fabs(EightValues[j]));

    }

    return CFL * h / MaxEihgtValue;
}

void Solver::Solve()
{
    double time = 0.0;
    double dt;
    int step = 0;
    while(time < T)
    {
        dt = compute_dt();
        HyperbolicOp.HyperbolicStep(dt);
        ++step;
    }

}