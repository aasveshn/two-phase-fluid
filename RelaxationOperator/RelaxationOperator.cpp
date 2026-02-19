#include "RelaxationOperator.h"

RelaxationOperator::RelaxationOperator(Mesh& msh, const Components& comp) : mesh(msh), 
                                                                            phases(comp){}

void RelaxationOperator::VelocityRelaxation(Cell& cell)
{

    double a2 = 1 - cell.W.a1;
    double UI = f_UI_x(cell.W);
    double VI = f_UI_y(cell.W);

    double e1_new = E(e_P_ro(cell.W.P1, cell.W.ro1, phases.p1), cell.W.u1, cell.W.v1) 
                    + UI*(UI-cell.W.u1) + VI*(VI-cell.W.v1) - UI*UI*0.5 - VI*VI*0.5;  
    double e2_new = E(e_P_ro(cell.W.P2, cell.W.ro2, phases.p2), cell.W.u2, cell.W.v2)
                    + UI*(UI-cell.W.u2) + VI*(VI-cell.W.v2) - UI*UI*0.5 - VI*VI*0.5; 
    
    cell.W.u1 = UI;
    cell.W.v1 = VI;
    cell.W.P1 = P_ro_e(cell.W.ro1, e1_new, phases.p1);
    cell.W.u2 = UI;
    cell.W.v2 = VI;
    cell.W.P2 = P_ro_e(cell.W.ro2, e2_new, phases.p2);
}

void RelaxationOperator::PressureRelaxation(Cell& cell)
{

    double a2 =  1 - cell.W.a1;
    double PI = f_PI(cell.W);
    const Phase& p1 = phases.p1;
    const Phase& p2 = phases.p2;
    StateW& W = cell.W;
    double dP = W.P1 - W.P2;

    

    if(dP == 0)
        return;

    double A0 = cell.W.a1*a2;
    double A1 =  A0 * (((p2.gamma + 1) * (W.P1 + p2.P0) + (p2.gamma-1)*(PI + p2.P0))/a2
                       + ((p1.gamma + 1) * (W.P2 + p1.P0) + (p1.gamma-1)*(PI + p1.P0))/ W.a1)*0.5;       
    double A2 =  (p2.gamma * (PI - p1.gamma * p1.P0)
                 - p1.gamma * (PI - p2.gamma*p2.P0) + (p2.gamma*p2.P0 - p1.gamma*p1.P0))*0.5;

    double A = 1.0;
    double B = A1 / A2;
    double C = -A0*dP/A2;
    double D = B*B - 4.0*A*C;
    double x1 = (-B + std::sqrt(D))*0.5;
    double x2 = (-B - std::sqrt(D))*0.5;
    double a1_0 = W.a1;

    if(std::fabs(x1) > std::fabs(x2))
        x2 = C / x1;
    else
        x1 = C / x2;

    double x;

    if(dP > 0)
        x1 > 0 ? x = x1 : x = x2;
    else 
    {
        if(a1_0 + x1 > 0 && a1_0 + x1 < 1)
            x = x1;
        else if (a1_0 + x2 > 0 && a1_0 + x2 < 1)
            x = x2;
        else
            std::cout<<"ERROR PRESSURE RELAXATION "<< x1 << " " << x2 << "\n"; 
   }
   

    
    W.a1 = a1_0 + x;

    double c1_1 = - ((p1.gamma - 1)*0.5*PI + p1.gamma*p1.P0);
    double c1_2 = - ((p2.gamma - 1)*0.5*PI + p2.gamma*p2.P0);
    double c2_1 = a1_0 * W.P1;
    double c2_2 = a2 * W.P2;
    double c3_1 = (p1.gamma + 1)*0.5;
    double c3_2 = (p2.gamma + 1)*0.5;

    W.P1 = (c1_1*(x) + c2_1) / (c3_1*(x) + a1_0);
    
    W.P2 = (c1_2*(-x) + c2_2) / (c3_2*(-x) + a2);

    W.ro1 = a1_0*W.ro1/(a1_0+x);
    W.ro2 = a2*W.ro2/(a2-x);

}



void RelaxationOperator::PressureTemperatureRelaxation(Cell& cell)
{
    StateW& W = cell.W;
    const Phase& p1 = phases.p1;
    const Phase& p2 = phases.p2;
    double a2 = 1 - W.a1;
    double PI = f_PI(W);

    double al_1 = p1.Cv*(p1.gamma-1)*W.a1*W.ro1;
    double al_2 = p2.Cv*(p2.gamma-1)*a2*W.ro2;
    
    double L = W.a1*((W.P1 + p1.gamma*p1.P0)/(p1.gamma-1)) + a2 * ((W.P2 + p2.gamma*p2.P0)/(p2.gamma-1));

    double a = (1.0 / (p1.gamma - 1) + (al_2/al_1) * (1.0/(p2.gamma-1)));
    double b = ((p2.P0 + p1.gamma*p1.P0)/(p1.gamma-1) - L + (al_2/al_1)*((p1.P0 + p2.gamma*p2.P0)/(p2.gamma-1) - L));
    double c = (p1.gamma*p1.P0*p2.P0)/(p1.gamma-1) - L*p2.P0 + (al_2/al_1)*((p2.gamma*p2.P0*p1.P0)/(p2.gamma-1) -L*p1.P0);
    double D = b*b - 4.0 * a * c;

    double x1, x2;

    if(D >= 0)
    {
        x1 = 0.5*(-b + std::sqrt(D))/a;
        x2 = 0.5*(-b - std::sqrt(D))/a;
    }
    else
    {
        std::cout<<"ERROR PRESS TEMP RELAX D < 0 ! \n";
        return;
    }

    double ak0 = p1.Cv*(p1.gamma - 1)*W.a1*W.ro1;
    double ak0_x1 = 1.0/(1.0 + (al_2/al_1) * ((x1 + p1.P0)/(x1 + p2.P0)));
    double ak0_x2 = 1.0/(1.0 + (al_2/al_1) * ((x2 + p1.P0)/(x2 + p2.P0)));

    double a1_x1 = ak0_x1;
    double a2_x1 = 1 - a1_x1;
    double a1_x2 = ak0_x2;
    double a2_x2 = 1 - a1_x2;

    double a1_new, P_eq;
    
    if((a1_x1 < 1 && a1_x1 > 0) && x1 > 0)
    {
        a1_new = a1_x1;
        P_eq = x1;
    }
    else if((a1_x2 < 1 && a1_x2 > 0) && x2 > 0)
    {
        a1_new = a1_x2;
        P_eq = x2;
    }
    else
    {
        std::cout<<"ERROR VOLUME FRACTION P&T RELAXATION \n";
        return;
    }

    W.ro1 = W.a1*W.ro1/a1_new;
    W.ro2 = a2*W.ro2/(1-a1_new);
    W.a1 = a1_new;
    W.P1 = P_eq;
    W.P2 = P_eq;
}

void RelaxationOperator::Relax()
{

    for(int i = 0; i < mesh.Cells.size(); ++i)
    {
        VelocityRelaxation(mesh.Cells[i]);
        PressureRelaxation(mesh.Cells[i]);
        if(mesh.Cells[i].is_Interface())
            PressureTemperatureRelaxation(mesh.Cells[i]);
    }
}