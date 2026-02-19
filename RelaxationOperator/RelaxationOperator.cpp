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


inline double RelaxationOperator::SolvePressure(double m1, double m2, double E, const Phase& p1, const Phase& p2)
{
    double A1 = m1*p1.Cv*(p1.gamma-1);
    double A2 = m2*p2.Cv*(p2.gamma-1);

    double a = m1*p1.Cv + m2*p2.Cv;
    double b = m1*p1.Cv*(p1.gamma*p1.P0 + p2.P0) + m2*p2.Cv*(p2.gamma*p2.P0 + p1.P0) -
               (E - m1*p1.q - m2*p2.q)*(A1+A2);
    double c = m1*p1.Cv*p1.gamma*p1.P0*p2.P0 + m2*p2.Cv*p2.gamma*p2.P0*p1.P0 - 
               (E - m1*p1.q - m2*p2.q)*(A1*p2.P0 + A2*p1.P0);
    
    double D = b*b - 4.0*a*c;
    
    double P1 = 0.5*(-b + std::sqrt(D))/a;
    double P2 = 0.5*(-b - std::sqrt(D))/a;

    double T1 = T_from_P(P1, m1, m2, phases);
    double T2 = T_from_P(P2, m1, m2, phases);

    if((P1 > 0 && T1 > 0) && (P2 > 0 && T2 > 0))
        std::cout<<"ERROR SOLVE PRESSURE GIBBS RELAXATION! \n";

    //?
    return (P1 > 0 && T1 > 0) ? P1 : P2;
}

inline double RelaxationOperator::Fdg(double m1, double m2, double E)
{
    double P = SolvePressure(m1, m2, E, phases.p1, phases.p2);
    double T = T_from_P(P, m1, m2, phases);

    double ro1 = ro_P_T(P,T,phases.p1);
    double ro2 = ro_P_T(P,T,phases.p2);

    double g1 = g_P_ro(P, ro1, phases.p1);
    double g2 = g_P_ro(P, ro2, phases.p2);

    return g1 - g2;   

}

inline double RelaxationOperator::Fdg(double P, double T)
{
    double ro1 = ro_P_T(P,T,phases.p1);
    double ro2 = ro_P_T(P,T,phases.p2);
   
    double g1 = g_P_ro(P, ro1, phases.p1);
    double g2 = g_P_ro(P, ro2, phases.p2);

    return g1 - g2;   

}

inline double RelaxationOperator::Fdg_der(double P, double T, const Phase& p1, const Phase& p2)
{
    double der_g1 = p1.gamma*p1.Cv - p1.qs - p1.Cv*p1.Cv*std::log(T) + p1.Cv*(p1.gamma-1)*std::log(P+p1.P0) - p1.Cv*p1.gamma;
    double der_g2 = p2.gamma*p2.Cv - p2.qs - p2.Cv*p2.Cv*std::log(T) + p2.Cv*(p2.gamma-1)*std::log(P+p2.P0) - p2.Cv*p2.gamma;
    return der_g1 - der_g2;
}

void RelaxationOperator::GibbsFreeEnergyRelaxation(Cell& cell)
{
    StateW& W = cell.W;
    const Phase& p1 = phases.p1;
    const Phase& p2 = phases.p2;
    double a2 = 1 - W.a1;

    double g1_0 = g_P_ro(W.P1, W.ro1, p1);
    double g2_0 = g_P_ro(W.P2, W.ro2, p2);

    double M = W.a1*W.ro1 + a2*W.ro2;
    double E = W.a1*W.ro1*e_P_ro(W.P1, W.ro1, p1) + a2*W.ro2*e_P_ro(W.P2, W.ro2, p2);

    double m1, m2;

    if(g1_0 - g2_0 > 0) // condensation
    {
        m1 = 1e-8;
        m2 = M - m2;

        double P = SolvePressure(m1, m2, E, p1, p2);
        double T = T_from_P(P, m1, m2, phases);

        double ro1 = ro_P_T(P, T, p1);
        double ro2 = ro_P_T(P, T, p2);

        double g1 = g_P_ro(P, ro1, p1);
        double g2 = g_P_ro(P, ro2, p2);

        double a1v = m1/ro1;
        double a2v = m2/ro2;

        double a1 = a1v/(a1v+a2v);

        if(g2 - g1 >= 0) // full condensation
        {
            W.a1 = a1;
            W.ro1 = ro1;
            W.P1 = P;
            W.ro2 = ro2;
            W.P2 = P;

            return;
        }
        else // partianl condensation
        {
            double m1_res, m2_res;

            double m1_min = 1e-8;
            double m1_max = W.a1*W.ro1;

            double m2_min = M - m1_min;
            double m2_max = M - m1_max;

            double f1 = Fdg(m1_min, m2_min, E);
            double f2 = Fdg(m1_max, m2_max, E);

            if(std::fabs(f1) < 1e-8)
            {
                m1_res = m1_min;
                m2_res = m2_min;
            }
            else if(std::fabs(f2) < 1e-8)
            {
                m1_res = m1_max;
                m2_res = m2_max;
            }
            else if(f1*f2 > 0)
            {
                if(std::fabs(f1) < std::fabs(f2))
                {
                    m1_res = m1_min;
                    m2_res = m2_min;
                }
                else
                {
                    m1_res = m1_max;
                    m2_res = m2_max;
                }
                
            }
            else
            {
                int maxIter = 1000;
                double tol = 1e-10;

                for(int i = 0; i < maxIter; ++i)
                {
                
                    double m1_mid = 0.5 * (m1_min + m1_max);
                    double m2_mid = M - m1_mid;
                    double fmid = Fdg(m1_mid, m2_mid, E);

                    if (std::fabs(fmid) < tol)
                        break;

                    if (f1 * fmid < 0.0)
                    {
                        m1_max  = m1_mid;
                        f2 = fmid;
                    }
                    else
                    {
                        m1_min = m1_mid;
                        f1 = fmid;
                    }
                }

                m1_res = 0.5 * (m1_min + m1_max);
                m2_res = M - m1_res;
            }

            double P_res = SolvePressure(m1_res, m2_res, E, p1, p2);
            double T_res = T_from_P(P_res, m1_res, m2_res, phases);

            double ro1_res = ro_P_T(P_res, T_res, p1);
            double ro2_res = ro_P_T(P_res, T_res, p2);

            double a1v = m1_res/ro1_res;
            double a2v = m2_res/ro2_res;
            double a1 = a1v/(a1v+a2v);

            W.a1 = a1;
            W.ro1 = ro1_res;
            W.ro2 = ro2_res;
            W.P1 = P_res;
            W.P2 = P_res;

            return;
        }

    }
    else if (fabs(g1_0 - g2_0) < 0.1) 
        return;
    else //evaporation
    {
        m1 = std::min(M-1e-8, (E-M*p2.q)/(p1.q-p2.q));
        m2 = M-m1;

        double P = SolvePressure(m1, m2, E, p1, p2);
        double T = T_from_P(P, m1, m2, phases);

        double ro1 = ro_P_T(P, T, p1);
        double ro2 = ro_P_T(P, T, p2);

        double g1 = g_P_ro(P, ro1, p1);
        double g2 = g_P_ro(P, ro2, p2);

        double a1v = m1/ro1;
        double a2v = m2/ro2;

        double a1 = a1v/(a1v+a2v);

        if(g1 - g2 <= 0) // full evaporation
        {
            W.a1 = a1;
            W.ro1 = ro1;
            W.P1 = P;
            W.ro2 = ro2;
            W.P2 = P;

            return;
        }
        else
        {
            double m1_res, m2_res;

            double m1_min = W.a1*W.ro1;
            double m1_max = m1;

            double m2_min = M - m1_min;
            double m2_max = M - m1_max;

            double f1 = Fdg(m1_min, m2_min, E);
            double f2 = Fdg(m1_max, m2_max, E);
            
             if(std::fabs(f1) < 1e-8)
            {
                m1_res = m1_min;
                m2_res = m2_min;
            }
            else if(std::fabs(f2) < 1e-8)
            {
                m1_res = m1_max;
                m2_res = m2_max;
            }
             if (f1 * f2 > 0.0)
            {
                   
                if(std::fabs(f1) < std::fabs(f2))
                {
                    m1_res = m1_min;
                    m2_res - m2_min;
                }
                else 
                {
                    m1_res = m1_max;
                    m2_res = m2_max;
                }
            }
            else
            {
                int maxIter = 1000;
                double tol = 1e-8;

                for(int i = 0; i < maxIter; ++i)
                {
                
                    double m1_mid = 0.5 * (m1_min + m1_max);
                    double m2_mid = M - m1_mid;
                    double fmid = Fdg(m1_mid, m2_mid, E);

                    if (std::fabs(fmid) < tol)
                        break;

                    if (f1 * fmid < 0.0)
                    {
                        m1_max  = m1_mid;
                        f2 = fmid;
                    }
                    else
                    {
                        m1_min = m1_mid;
                        f1 = fmid;
                    }
                }

                m1_res = 0.5 * (m1_min + m1_max);
                m2_res = M - m1_res;
            }

            double P_res = SolvePressure(m1_res, m2_res, E, p1, p2);
            double T_res = T_from_P(P_res, m1_res, m2_res, phases);

            double ro1_res = ro_P_T(P_res, T_res, p1);
            double ro2_res = ro_P_T(P_res, T_res, p2);

            double a1v = m1_res/ro1_res;
            double a2v = m2_res/ro2_res;
            double a1 = a1v/(a1v+a2v);

            W.a1 = a1;
            W.ro1 = ro1_res;
            W.ro2 = ro2_res;
            W.P1 = P_res;
            W.P2 = P_res;

            return;
        }
    }

}


double RelaxationOperator::SolveTemperature(double P, double T0, double tol, int maxIter)
{

    double T = T0;
    for(int i = 200; i < 1000; ++i)
    {
        if(Fdg(P, i)*Fdg(P, i + 10) < 0)
        {
            T = i+5;
            break;
        }
    }
    
    
    for (int iter = 0; iter < maxIter; ++iter) 
    {
        double f = Fdg(P, T);
        double df = Fdg_der(P, T, phases.p1, phases.p2);
        double delta = f / df;
        T -= delta;
        if (std::fabs(delta) < tol) {


            return T;
        }
    }
    std::cout << T<<" "<<T0 <<" "<<P<<"\n";
    throw std::runtime_error("method did not converge within the maximum number of iterations.");
}


void RelaxationOperator::Relax()
{

    for(int i = 0; i < mesh.Cells.size(); ++i)
    {
        VelocityRelaxation(mesh.Cells[i]);
        PressureRelaxation(mesh.Cells[i]);
        if(mesh.Cells[i].is_Interface())
            PressureTemperatureRelaxation(mesh.Cells[i]);
       if(mesh.Cells[i].is_Interface())
        {
           double T = T_ro_P(mesh.Cells[i].W.ro1, mesh.Cells[i].W.P1, phases.p1);
            //точность?
            double Tsat = SolveTemperature(mesh.Cells[i].W.P1, T, 1e-2, 1000);
            
            if(T > Tsat)
            {
               GibbsFreeEnergyRelaxation(mesh.Cells[i]);
            } 
        }
    }
}