#include "RiemanSolver.h"

RiemanSolver::RiemanSolver(unsigned int N, unsigned int k): Flux_X(k), Flux_Y(N-k){};

void RiemanSolver::HLLC(const Mesh& mesh, const std::vector<StateW>& WL_in, const std::vector<StateW>& WR_in,  
                        unsigned int H_idx, bool is_X_dir, const Components& phases)
{
    
    unsigned int N, interval_L, interval_R;
    unsigned int face_dirs[2];
    std::vector<StateU>& flux = is_X_dir ? Flux_X : Flux_Y;

    if(is_X_dir)
    {
        interval_L = 0;
        interval_R = H_idx -1;
        N = H_idx;
        face_dirs[0] = 0;
        face_dirs[1] = 1;
    }
    else
    {
        interval_L = H_idx;
        interval_R = mesh.Faces.size()-1;
        N = mesh.Faces.size() - H_idx;
        face_dirs[0] = 2;
        face_dirs[1] = 3;
    }

    int idx_L, idx_R;
    
    StateU UL_C;
    StateU UR_C;

    for(int i = interval_L; i < N; ++i)
    {

        
        idx_L = mesh.Faces[i].Left_ID;
        idx_R = mesh.Faces[i].Right_ID;

        if(idx_L == -1 || idx_R == -1)
        {    
            flux[i].zero();
            continue;
        }

        const StateW& WL = WR_in[idx_L];
        const StateW& WR = WL_in[idx_R];

        double a2L = 1.0 - WL.a1;
        double a2R = 1.0 - WR.a1;
        
        StateU UL(WL, phases);
        StateU UR(WR, phases);

        double C1L = f_C(WL.P1, WL.ro1, phases.p1);
        double C2L = f_C(WL.P2, WL.ro2, phases.p2);

        double C1R = f_C(WR.P1, WR.ro1, phases.p1);
        double C2R = f_C(WR.P2, WR.ro2, phases.p2);
        

        double sL = std::min({WL.u1 - C1L, WL.u2 - C2L, WR.u1 - C1R, WR.u2 - C2R});
        double sR = std::max({WL.u1 + C2L, WL.u2 + C2L, WR.u1 + C1R, WR.u2 + C2R});


        double roL = WL.a1*WL.ro1 + a2L*WL.ro2;
        double PL = WL.a1*WL.P1 + a2L*WL.P2;
        double uL = (WL.a1*WL.ro1*WL.u1 + a2L*WL.ro2*WL.u2)/roL;

        double roR = WR.a1*WR.ro1 + a2R*WR.ro2;
        double PR = WR.a1*WR.P1 + a2R*WR.P2;
        double uR = (WR.a1*WR.ro1*WR.u1 + a2R*WR.ro2*WR.u2)/roR;

        double s_C = (PR - PL + roL*uL*(sL-uL) - roR*uR*(sR-uR))/(roL*(sL-uL) - roR*(sR-uR));
        double d_sL = sL - s_C;
        double d_sR = sR - s_C;

       

        UL_C[0] = WL.a1*WL.ro1*((sL-WL.u1)/d_sL);
        UL_C[1] = UL_C[0]*s_C;
        UL_C[2] = UL_C[0]*(E(e_P_ro(WL.P1, WL.ro1, phases.p1), WL.u1)+
                  (s_C-WL.u1)*(s_C + WL.P1/(WL.ro1*(sL-WL.u1))));
        UL_C[3] = a2L*WL.ro2*((sL-WL.u2)/d_sL);
        UL_C[4] = UL_C[0]*s_C;
        UL_C[5] = UL_C[0]*(E(e_P_ro(WL.P2, WL.ro2, phases.p2), WL.u2)+
                  (s_C-WL.u2)*(s_C + WL.P2/(WL.ro2*(sL-WL.u2))));

        UR_C[0] = WR.a1*WR.ro1*((sR-WR.u1)/d_sR);
        UR_C[1] = UR_C[0]*s_C;
        UR_C[2] = UR_C[0]*(E(e_P_ro(WR.P1, WR.ro1, phases.p1), WR.u1)+
                  (s_C-WR.u1)*(s_C + WR.P1/(WR.ro1*(sR-WR.u1))));
        UR_C[3] = a2R*WR.ro2*((sR-WR.u2)/d_sR);
        UR_C[4] = UR_C[0]*s_C;
        UR_C[5] = UR_C[0]*(E(e_P_ro(WR.P2, WR.ro2, phases.p2), WR.u2)+
                  (s_C-WR.u2)*(s_C + WR.P2/(WR.ro2*(sR-WR.u2))));

        if(sL >= 0)
        {
            flux[i][0] = UL[1];
            flux[i][1] = UL[1]*WL.u1 + WL.a1*WL.P1;
            flux[i][2] = WL.a1*(WL.ro1*(UL[2]/UL[0]) + WL.P1)*WL.u1;
            flux[i][3] = UL[4];
            flux[i][4] = UL[4]*WL.u2 + a2L*WL.P2;
            flux[i][5] = a2L*(WL.ro2*(UL[5]/UL[3]) + WL.P2)*WL.u2;
        }
        else if(sL <= 0 && s_C >= 0)
        {
            flux[i][0] = UL[1] + sL*(UL_C[0] - UL[0]);
            flux[i][1] = UL[1] * WL.u1 + WL.a1*WL.P1 + sL*(UL_C[1] - UL[1]);
            flux[i][2] = WL.a1 * (WL.ro1*(UL[2]/UL[0]) + WL.P1) * WL.u1 + sL*(UL_C[2] - UL[2]);
            flux[i][3] = UL[4] + sL*(UL_C[3] - UL[3]);
            flux[i][4] = UL[4] * WL.u2 + a2L*WL.P2 + sL*(UL_C[4] - UL[4]);
            flux[i][5] = a2L * (WL.ro2*(UL[5]/UL[3]) + WL.P2) * WL.u2 + sL*(UL_C[5] - UL[5]);
        }
        else if(s_C <= 0 && sR >= 0)
        {
            flux[i][0] = UR[1] + sR*(UR_C[0] - UR[0]);
            flux[i][1] = UR[1] * WR.u1 + WR.a1*WR.P1 + sR*(UR_C[1] - UR[1]);
            flux[i][2] = WR.a1 * (WR.ro1*(UR[2]/UR[0]) + WR.P1) * WR.u1 + sR*(UR_C[2] - UR[2]);
            flux[i][3] = UR[4] + sR*(UR_C[3] - UR[3]);
            flux[i][4] = UR[4] * WR.u2 + a2R*WR.P2 + sR*(UR_C[4] - UR[4]);
            flux[i][5] = a2R * (WR.ro2*(UR[5]/UR[3]) + WR.P2) * WR.u2 + sR*(UR_C[5] - UR[5]);
        }
        else if(sR <= 0)
        {
            flux[i][0] = UR[1];
            flux[i][1] = UR[1]*WR.u1 + WR.a1*WR.P1;
            flux[i][2] = WR.a1*(WR.ro1*(UR[2]/UR[0]) + WR.P1)*WR.u1;
            flux[i][3] = UR[4];
            flux[i][4] = UR[4]*WR.u2 + a2R*WR.P2;
            flux[i][5] = a2R*(WR.ro2*(UR[5]/UR[3]) + WR.P2)*WR.u2;
        }
        else
        {
            std::cout<<"HLLC ERROR \n";
        }
    }
}

