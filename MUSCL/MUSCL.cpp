#include "MUSCL.h"

MUSCL::MUSCL(int N) : B(betta), W_L(N), W_R(N){};

StateW MUSCL::slopLimiter(const StateW& dL, const StateW& dR)
{
    StateW result;

    for(int i = 0; i < 9; ++i)
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

StateW MUSCL::A_prod_W(const StateW& W, bool is_X_dir, const StateW& WL, const StateW& WR, const Components& phases)
{
    StateW result;
    StateW dW = WR - WL;

    double a1 = W.a1;
    double a2 = 1.0 - a1;

    double u1_n = is_X_dir ? W.u1 : W.v1;
    double u2_n = is_X_dir ? W.u2 : W.v2;
    
    double UI = is_X_dir ? f_UI_x(W) : f_UI_y(W); 
    double PI = f_PI(W); 

    double CI1_sq = f_C2(PI, W.ro1, phases.p1); // c_I^2
    double CI2_sq = f_C2(PI, W.ro2, phases.p2);
    double C1_sq  = f_C2(W.P1, W.ro1, phases.p1); // c^2
    double C2_sq  = f_C2(W.P2, W.ro2, phases.p2);


    result[0] = UI * dW[0];

    if (is_X_dir) {
       
        result[1] = (W.ro1 * (u1_n - UI) / a1) * dW[0] + u1_n * dW[1] + W.ro1 * dW[2];
        result[2] = ((W.P1 - PI) / (a1 * W.ro1)) * dW[0] + u1_n * dW[2] + (1.0 / W.ro1) * dW[4];
        result[3] = u1_n * dW[3];
        result[4] = (W.ro1 * CI1_sq * (u1_n - UI) / a1) * dW[0] + (W.ro1 * C1_sq) * dW[2] + u1_n * dW[4];
        result[5] = (-W.ro2 * (u2_n - UI) / a2) * dW[0] + u2_n * dW[5] + W.ro2 * dW[6];
        result[6] = (-(W.P2 - PI) / (a2 * W.ro2)) * dW[0] + u2_n * dW[6] + (1.0 / W.ro2) * dW[8];
        result[7] = u2_n * dW[7];
        result[8] = (-W.ro2 * CI2_sq * (u2_n - UI) / a2) * dW[0] + (W.ro2 * C2_sq) * dW[6] + u2_n * dW[8];
    } 
    else {

        result[1] = (W.ro1 * (u1_n - UI) / a1) * dW[0] + u1_n * dW[1] + W.ro1 * dW[3]; 
        result[2] = u1_n * dW[2];
        result[3] = ((W.P1 - PI) / (a1 * W.ro1)) * dW[0] + u1_n * dW[3] + (1.0 / W.ro1) * dW[4];
        result[4] = (W.ro1 * CI1_sq * (u1_n - UI) / a1) * dW[0] + (W.ro1 * C1_sq) * dW[3] + u1_n * dW[4];
        result[5] = (-W.ro2 * (u2_n - UI) / a2) * dW[0] + u2_n * dW[5] + W.ro2 * dW[7];
        result[6] = u2_n * dW[6];
        result[7] = (-(W.P2 - PI) / (a2 * W.ro2)) * dW[0] + u2_n * dW[7] + (1.0 / W.ro2) * dW[8];
        result[8] = (-W.ro2 * CI2_sq * (u2_n - UI) / a2) * dW[0] + (W.ro2 * C2_sq) * dW[7] + u2_n * dW[8];
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
            d_L = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            d_R = mesh.Cells[R_id].W - mesh.Cells[i].W;
        }
        else if(mesh.Faces[mesh.Cells[i].faces_ID[face_dirs[1]]].type != 0)
        {
            L_id = mesh.Faces[mesh.Cells[i].faces_ID[face_dirs[0]]].Left_ID;
            d_L = mesh.Cells[i].W - mesh.Cells[L_id].W;
            d_R = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
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
               
        AWprod = A_prod_W(mesh.Cells[i].W, is_X_dir, W_L[i], W_R[i], phases);

        W_L[i] = W_L[i] - (0.5*dt/h) * AWprod;
        W_R[i] = W_R[i] - (0.5*dt/h) * AWprod;
        
    }

}

