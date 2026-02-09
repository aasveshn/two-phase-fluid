#include "HyperbolicOperator.h"

HyperbolicOperator::HyperbolicOperator(Mesh& msh, const Components& comp) :phases(comp),
                                                                           muscl(msh.Cells.size()),
                                                                           rieamn(msh.Faces.size(), msh.VertToHoriz),
                                                                           mesh(msh){}


void HyperbolicOperator::GodunovStep(bool is_X_dir, double dt)
{
    int N = mesh.Cells.size();
    unsigned int face_dirs[2];
    double h;
    unsigned int k;
    if(is_X_dir)
    {
        face_dirs[0] = 0;
        face_dirs[1] = 1;
        h = dx;
        k = 0;
    }
    else
    {
        face_dirs[0] = 2;
        face_dirs[1] = 3;
        h = dy;
        k = mesh.VertToHoriz;
    }
    std::vector<StateU> flux = is_X_dir ? rieamn.getXFlux() : rieamn.getYFlux();
    StateW& W = mesh.Cells[0].W;
    unsigned int face_L, face_R;

    for(int i = 0; i < N; ++i)
    {

        W = mesh.Cells[i].W;
        StateU U_prev(W, phases);
        StateU U_new;

        double a2 = 1.0 - W.a1;
        double UI = f_UI(W);
        double PI = f_PI(W); 

        double da_1 = muscl.getW_R()[i].a1 - muscl.getW_L()[i].a1;

        face_L = mesh.Cells[i].faces_ID[face_dirs[0]] - k;
        face_R = mesh.Cells[i].faces_ID[face_dirs[1]] - k;

        //StateU& flux_L = 
        //StateU& flux_R;

        U_new[0] = U_prev[0] - (dt/h)*(flux[face_R][0] - flux[face_L][0]);
        U_new[1] = U_prev[1] - (dt/h)*((flux[face_R][1] - flux[face_L][1]) +  da_1*PI);
        U_new[2] = U_prev[2] - (dt/h)*((flux[face_R][2] - flux[face_L][2]) +  da_1*PI*UI);
        U_new[3] = U_prev[3] - (dt/h)*(flux[face_R][3] - flux[face_L][3]);
        U_new[4] = U_prev[4] - (dt/h)*((flux[face_R][4] - flux[face_L][4]) -  da_1*PI);
        U_new[5] = U_prev[5] - (dt/h)*((flux[face_R][5] - flux[face_L][5]) -  da_1*PI*UI);

        double a1_new = mesh.Cells[i].W.a1 - UI*(dt/h)*(da_1);
        double a2_new = 1 - a1_new;

        double E1_new = U_new[2]/U_new[0];
        double E2_new = U_new[5]/U_new[3];
        

        mesh.Cells[i].W.a1 = a1_new;
        mesh.Cells[i].W.ro1 = U_new[0]/a1_new;
        mesh.Cells[i].W.u1 = U_new[1]/U_new[0];
        double e1_new = E1_new - mesh.Cells[i].W.u1*mesh.Cells[i].W.u1*0.5;
        mesh.Cells[i].W.P1 = P_ro_e(mesh.Cells[i].W.ro1, e1_new, phases.p1);

        mesh.Cells[i].W.ro2 = U_new[3]/a2_new;
        mesh.Cells[i].W.u2 = U_new[4]/U_new[3];
        double e2_new = E2_new - mesh.Cells[i].W.u2*mesh.Cells[i].W.u2*0.5;
        mesh.Cells[i].W.P2 = P_ro_e(mesh.Cells[i].W.ro2, e2_new, phases.p2);
    }

}
                                                                
void HyperbolicOperator::HyperbolicStep(bool is_X_dir, double dt)
{
    //X step
    muscl.MUSCL_Operator(mesh, true, dt, phases);
    rieamn.HLLC(mesh, muscl.getW_L(), muscl.getW_R(), mesh.VertToHoriz, true, phases);
    GodunovStep(true, dt);
    //Y strep
    muscl.MUSCL_Operator(mesh, false, dt, phases);
    rieamn.HLLC(mesh, muscl.getW_L(), muscl.getW_R(), mesh.VertToHoriz, false, phases);
    GodunovStep(false, dt);

}                                                                           