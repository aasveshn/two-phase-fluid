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
                           HyperbolicOp(msh, phases),
                           RelaxationOp(msh, phases){}

double Solver::compute_dt()
{
    double MaxValue = -1.0; 
    double h = std::max(dx, dy);

    for (int i = 0; i < mesh.Cells.size(); ++i)
    {
        const StateW& W = mesh.Cells[i].W;
        
        double C1 = f_C(W.P1, W.ro1, phases.p1);
        double C2 = f_C(W.P2, W.ro2, phases.p2);
        
        double UIx = f_UI_x(W);
        double UIy = f_UI_y(W);

     
        double cellMax = std::max({
            
            std::abs(UIx),              
            std::abs(W.u1 - C1),       
            std::abs(W.u1),             
            std::abs(W.u1),             
            std::abs(W.u1 + C1),       
            std::abs(W.u2 - C2),      
            std::abs(W.u2),            
            std::abs(W.u2),             
            std::abs(W.u2 + C2),        

           
            std::abs(UIy),              
            std::abs(W.v1 - C1),       
            std::abs(W.v1),            
            std::abs(W.v1),            
            std::abs(W.v1 + C1),        
            std::abs(W.v2 - C2),       
            std::abs(W.v2),             
            std::abs(W.v2),             
            std::abs(W.v2 + C2)         
        });

        if (cellMax > MaxValue) MaxValue = cellMax;
    }

    return CFL * h / MaxValue;
}

void Solver::Solve()
{
    //Делегировать запись в файлы в отедльный метод
    double time = 0.0;
    double dt = 1e-10;
    int step = 0;

    std::ofstream f_a1("output/a1.txt");
    std::ofstream f_ro1("output/ro1.txt");
    std::ofstream f_u1("output/u1.txt");
    std::ofstream f_v1("output/v1.txt");
    std::ofstream f_P1("output/P1.txt");
    std::ofstream f_ro2("output/ro2.txt");
    std::ofstream f_u2("output/u2.txt");
    std::ofstream f_v2("output/v2.txt");
    std::ofstream f_P2("output/P2.txt");

     if (!f_a1.is_open()) {
        std::cerr << "Error opening output files!" << std::endl;
        return;
    }

   
   //1500   
    while(time < T)
    {
        HyperbolicOp.HyperbolicStep(dt);
        time += dt;
        dt = compute_dt();
        std::cout<<time<<"\n";
        ++step;
        if(step == 50)
            break;
    }


     for (const auto& cell : mesh.Cells) {
        f_a1  << cell.W.a1  << "\n";
        f_ro1 << cell.W.ro1 << "\n";
        f_u1  << cell.W.u1  << "\n";
        f_v1  << cell.W.v1  << "\n";
        f_P1  << cell.W.P1  << "\n";
        f_ro2 << cell.W.ro2 << "\n";
        f_u2  << cell.W.u2  << "\n";
        f_v2  << cell.W.v2  << "\n";
        f_P2  << cell.W.P2  << "\n";
    }

    f_a1.close(); f_ro1.close(); f_u1.close(); f_v1.close(); f_P1.close();
    f_ro2.close(); f_u2.close(); f_v2.close(); f_P2.close();
    
}