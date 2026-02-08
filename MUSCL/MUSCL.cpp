#include "MUSCL.h"

MUSCL::MUSCL(int N) : B(betta), W_L(N), W_R(N){};
//MUSCL::MUSCL(int N) : B(betta), W_L(N), W_R(N), W_ex_L(N), W_ex_R(N){};



StateW MUSCL::slopLimiter(const StateW& dL, const StateW& dR)
{
    StateW result;

    for(int i = 0; i < 7; ++i)
    {
        if(dR[i] > 0)
            result[i] = std::max({0.0, 
                                  std::min(B*dL[i], dR[i]),
                                  std::min(dL[i], B*dR[i])});
        else
            result[i] = std::min({0.0, 
                                  std::max(B*dL[i], dR[i]),
                                  std::max(dL[i], B*dR[i])});    
    }

    return result;
}

StateW MUSCL::A_prod_W(const StateW& W, const StateW& WL, const StateW& WR, const Components& phases)
{
    StateW result;
    std::vector<std::vector<double>> A(7, std::vector<double>(7,0));
    double a2 = 1.0 - W.a1;
    double UI = f_UI(W);
    double PI = f_PI(W);
    double CI1 = f_C(PI, W.ro1, phases.p1);
    double CI2 = f_C(PI, W.ro2, phases.p2);
    double C1 = f_C(W.P1, W.ro1, phases.p1);
    double C2 = f_C(W.P2, W.ro2, phases.p2);

    A[0][0] = UI;

    A[1][0] = W.ro1*(W.u1 - UI)/W.a1;
    A[1][1] = W.u1;
    A[1][2] = W.ro1;

    A[2][0] = (W.P1 - PI)/(W.a1*W.ro1);
    A[2][2] = W.u1;
    A[2][3] = 1.0 / W.ro1;

    A[3][0] = W.ro1*CI1*CI1*(W.u1 - UI) / W.a1;
    A[3][2] = W.ro1 * C1*C1;
    A[3][3] = W.u1;

    A[4][0] = - W.ro2*(W.u2 - UI)/a2;
    A[4][4] = W.u2;
    A[4][5] = W.ro2;

    A[5][0] = - (W.P2 - PI)/(a2*W.ro2);
    A[5][5] = W.u2;
    A[5][6] = 1.0 / W.ro2;

    A[6][0] = - (W.ro2*CI2*CI2)*(W.u2 - UI)/a2;
    A[6][5] = W.ro2*C2*C2;
    A[6][6] = W.u2;

    for(int i = 0; i < 7; ++i)
    {
        for(int j = 0; j < 0; ++j)
            result[i] += A[i][j]*(WR[j]-WL[j]);

    }

    return result;
}

void MUSCL::MUSCL_Operator(const Mesh& mesh, bool is_X_dir, double dt, const Components& phases)
{
    double h;
    unsigned int face_dirs[2];
    if(is_X_dir)
    {
        face_dirs[0] = 0;
        face_dirs[1] = 1;
        h = dx;
    }
    else
    {
        face_dirs[0] = 2;
        face_dirs[1] = 3;
        h = dy;
    }
    
    int N = W_L.size();

    
    unsigned R_id, L_id;
    StateW d_L, d_R, slope, AWprod;
    for(int i = 0; i < N; ++i)
    {
        
        if(mesh.Faces[mesh.Cells[i].faces_ID[face_dirs[0]]].type != 0)
        {
            R_id = mesh.Faces[mesh.Cells[i].faces_ID[face_dirs[1]]].Right_ID;
            d_L = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            d_R = mesh.Cells[R_id].W - mesh.Cells[i].W;
        }
        if(mesh.Faces[mesh.Cells[i].faces_ID[face_dirs[1]]].type != 0)
        {
            L_id = mesh.Faces[mesh.Cells[i].faces_ID[face_dirs[0]]].Left_ID;
            d_L = mesh.Cells[i].W - mesh.Cells[L_id].W;
            d_R = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        }
        else
        {
            R_id = mesh.Faces[mesh.Cells[i].faces_ID[face_dirs[1]]].Right_ID;
            L_id = mesh.Faces[mesh.Cells[i].faces_ID[face_dirs[0]]].Left_ID;
            d_L = mesh.Cells[i].W - mesh.Cells[L_id].W;
            d_R = mesh.Cells[R_id].W - mesh.Cells[i].W;
        }

        slope = slopLimiter(d_L, d_R);

        W_L[i] = mesh.Cells[i].W - 0.5 * slope;
        W_R[i] = mesh.Cells[i].W + 0.5 * slope;
                
        AWprod = A_prod_W(mesh.Cells[i].W, W_L[i], W_R[i], phases);

        W_L[i] = W_L[i] - (0.5*dt/h) * AWprod;
        W_R[i] = W_R[i] - (0.5*dt/h) * AWprod;
    }

}

