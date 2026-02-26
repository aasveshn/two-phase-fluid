#include "RiemanSolver.h"

RiemanSolver::RiemanSolver(unsigned int N, unsigned int k): Flux_X(k+1), Flux_Y(N-k){};

void RiemanSolver::HLLC(const Mesh& mesh, const std::vector<StateW>& WL_in, const std::vector<StateW>& WR_in,  
                        unsigned int H_idx, bool is_X_dir, const Components& phases)
{
    unsigned int interval_L, N;
    int k;
    std::vector<StateU>& flux = is_X_dir ? Flux_X : Flux_Y;

    if(is_X_dir) {
        interval_L = 0;
        N = H_idx+1; 
        k = 0;
    } else {
        interval_L = H_idx;
        N = mesh.Faces.size(); 
        k = H_idx;
    }

    for(unsigned int i = interval_L; i < N; ++i)
    {

        int idx_L = mesh.Faces[i].Left_ID;
        int idx_R = mesh.Faces[i].Right_ID;

        StateW WL;
        StateW WR;
 
        if(idx_L == -1)
        {
            idx_L = idx_R;
            WL = WL_in[idx_L];
            WR = WL_in[idx_R];
            
            
        } 
        else if(idx_R == -1) 
        {
            idx_R = idx_L;
            WR = WR_in[idx_R];
            WL = WR_in[idx_L];
        }
        else
        {
            WL = WR_in[idx_L];
            WR = WL_in[idx_R];
        }
        

        

        if(mesh.Faces[i].type == 1)
        {
            if(mesh.Faces[i].Left_ID == -1)
            {
                
                if(is_X_dir)
                {
                    WL.u1 = -WL.u1;
                    WL.u2 = -WL.u2;
                }
                else
                {
                    WL.v1 = -WL.v1;
                    WL.v2 = -WL.v2;
                }
            }
            else
            {
                
                if(is_X_dir)
                {
                    WR.u1 = -WR.u1;
                    WR.u2 = -WR.u2;
                }
                else
                {
                    WR.v1 = -WR.v1;
                    WR.v2 = -WR.v2;
                }
            }

        }
        if(mesh.Faces[i].type == 2)
        {
            WL = StateW(a1L0, ro1L0, u1L0, v1L0, P1L0, ro2L0, u2L0, v2L0, P2L0);
        }

        double a1L = WL.a1; double a2L = 1.0 - a1L;
        double a1R = WR.a1; double a2R = 1.0 - a1R;
        
        StateU UL(WL, phases); 
        StateU UR(WR, phases);

        double C1L = f_C(WL.P1, WL.ro1, phases.p1);
        double C2L = f_C(WL.P2, WL.ro2, phases.p2);
        double C1R = f_C(WR.P1, WR.ro1, phases.p1);
        double C2R = f_C(WR.P2, WR.ro2, phases.p2);
        
        double roL = a1L*WL.ro1 + a2L*WL.ro2;
        double PL  = a1L*WL.P1  + a2L*WL.P2;
        double roR = a1R*WR.ro1 + a2R*WR.ro2;
        double PR  = a1R*WR.P1  + a2R*WR.P2;
        
        double sL, sR, uL, uR;
        StateU UL_C, UR_C;

        if(is_X_dir)
        {
            
            sL = std::min({WL.u1 - C1L, WL.u2 - C2L, WR.u1 - C1R, WR.u2 - C2R});
            sR = std::max({WL.u1 + C1L, WL.u2 + C2L, WR.u1 + C1R, WR.u2 + C2R});
            uL = (a1L*WL.ro1*WL.u1 + a2L*WL.ro2*WL.u2)/roL;
            uR = (a1R*WR.ro1*WR.u1 + a2R*WR.ro2*WR.u2)/roR;

            double s_C = (PR - PL + roL*uL*(sL-uL) - roR*uR*(sR-uR))/(roL*(sL-uL) - roR*(sR-uR));
            
           
            UL_C[0] = a1L*WL.ro1*((sL-WL.u1)/(sL-s_C));
            UL_C[1] = UL_C[0]*s_C;   
            UL_C[2] = UL_C[0]*WL.v1; 
            UL_C[3] = UL_C[0]*(E(e_P_ro(WL.P1, WL.ro1, phases.p1), WL.u1, WL.v1) + (s_C-WL.u1)*(s_C + WL.P1/(WL.ro1*(sL-WL.u1))));
            
            
            UL_C[4] = a2L*WL.ro2*((sL-WL.u2)/(sL-s_C));
            UL_C[5] = UL_C[4]*s_C;   
            UL_C[6] = UL_C[4]*WL.v2; 
            UL_C[7] = UL_C[4]*(E(e_P_ro(WL.P2, WL.ro2, phases.p2), WL.u2, WL.v2) + (s_C-WL.u2)*(s_C + WL.P2/(WL.ro2*(sL-WL.u2))));

            
            UR_C[0] = a1R*WR.ro1*((sR-WR.u1)/(sR-s_C));
            UR_C[1] = UR_C[0]*s_C; UR_C[2] = UR_C[0]*WR.v1;
            UR_C[3] = UR_C[0]*(E(e_P_ro(WR.P1, WR.ro1, phases.p1), WR.u1, WR.v1) + (s_C-WR.u1)*(s_C + WR.P1/(WR.ro1*(sR-WR.u1))));
            UR_C[4] = a2R*WR.ro2*((sR-WR.u2)/(sR-s_C));
            UR_C[5] = UR_C[4]*s_C; UR_C[6] = UR_C[4]*WR.v2;
            UR_C[7] = UR_C[4]*(E(e_P_ro(WR.P2, WR.ro2, phases.p2), WR.u2, WR.v2) + (s_C-WR.u2)*(s_C + WR.P2/(WR.ro2*(sR-WR.u2))));

            
            if(sL >= 0) {
                flux[i-k][0] = UL[1]; 
                flux[i-k][1] = UL[1]*WL.u1 + a1L*WL.P1; 
                flux[i-k][2] = UL[1]*WL.v1; 
                flux[i-k][3] = (UL[3] + a1L*WL.P1)*WL.u1;
                flux[i-k][4] = UL[5];
                flux[i-k][5] = UL[5]*WL.u2 + a2L*WL.P2; 
                flux[i-k][6] = UL[5]*WL.v2; 
                flux[i-k][7] = (UL[7] + a2L*WL.P2)*WL.u2;
            } else if(sL <= 0 && s_C >= 0) {
                for(int j=0; j<8; ++j) flux[i-k][j] = get_F_X(UL, WL, a1L, a2L)[j] + sL*(UL_C[j] - UL[j]);
            } else if(s_C <= 0 && sR >= 0) {
                for(int j=0; j<8; ++j) flux[i-k][j] = get_F_X(UR, WR, a1R, a2R)[j] + sR*(UR_C[j] - UR[j]);
            } else {
                flux[i-k][0] = UR[1]; 
                flux[i-k][1] = UR[1]*WR.u1 + a1R*WR.P1;
                flux[i-k][2] = UR[1]*WR.v1; 
                flux[i-k][3] = (UR[3] + a1R*WR.P1)*WR.u1;
                flux[i-k][4] = UR[5];
                flux[i-k][5] = UR[5]*WR.u2 + a2R*WR.P2; 
                flux[i-k][6] = UR[5]*WR.v2;
                flux[i-k][7] = (UR[7] + a2R*WR.P2)*WR.u2;
            }
        }
        else
        {
           
            sL = std::min({WL.v1 - C1L, WL.v2 - C2L, WR.v1 - C1R, WR.v2 - C2R});
            sR = std::max({WL.v1 + C1L, WL.v2 + C2L, WR.v1 + C1R, WR.v2 + C2R});
            uL = (a1L*WL.ro1*WL.v1 + a2L*WL.ro2*WL.v2)/roL;
            uR = (a1R*WR.ro1*WR.v1 + a2R*WR.ro2*WR.v2)/roR;

            double s_C = (PR - PL + roL*uL*(sL-uL) - roR*uR*(sR-uR))/(roL*(sL-uL) - roR*(sR-uR));

            
            UL_C[0] = a1L*WL.ro1*((sL-WL.v1)/(sL-s_C));
            UL_C[1] = UL_C[0]*WL.u1; 
            UL_C[2] = UL_C[0]*s_C;   
            UL_C[3] = UL_C[0]*(E(e_P_ro(WL.P1, WL.ro1, phases.p1), WL.u1, WL.v1) + (s_C-WL.v1)*(s_C + WL.P1/(WL.ro1*(sL-WL.v1))));

           
            UL_C[4] = a2L*WL.ro2*((sL-WL.v2)/(sL-s_C));
            UL_C[5] = UL_C[4]*WL.u2; 
            UL_C[6] = UL_C[4]*s_C;   
            UL_C[7] = UL_C[4]*(E(e_P_ro(WL.P2, WL.ro2, phases.p2), WL.u2, WL.v2) + (s_C-WL.v2)*(s_C + WL.P2/(WL.ro2*(sL-WL.v2))));

            
            UR_C[0] = a1R*WR.ro1*((sR-WR.v1)/(sR-s_C));
            UR_C[1] = UR_C[0]*WR.u1; UR_C[2] = UR_C[0]*s_C;
            UR_C[3] = UR_C[0]*(E(e_P_ro(WR.P1, WR.ro1, phases.p1), WR.u1, WR.v1) + (s_C-WR.v1)*(s_C + WR.P1/(WR.ro1*(sR-WR.v1))));
            UR_C[4] = a2R*WR.ro2*((sR-WR.v2)/(sR-s_C));
            UR_C[5] = UR_C[4]*WR.u2; UR_C[6] = UR_C[4]*s_C;
            UR_C[7] = UR_C[4]*(E(e_P_ro(WR.P2, WR.ro2, phases.p2), WR.u2, WR.v2) + (s_C-WR.v2)*(s_C + WR.P2/(WR.ro2*(sR-WR.v2))));

            if(sL >= 0) {
                flux[i-k][0] = UL[2];
                flux[i-k][1] = UL[2]*WL.u1; 
                flux[i-k][2] = UL[2]*WL.v1 + a1L*WL.P1; 
                flux[i-k][3] = (UL[3] + a1L*WL.P1)*WL.v1;
                flux[i-k][4] = UL[6]; 
                flux[i-k][5] = UL[6]*WL.u2; 
                flux[i-k][6] = UL[6]*WL.v2 + a2L*WL.P2; 
                flux[i-k][7] = (UL[7] + a2L*WL.P2)*WL.v2;
            } else if(sL <= 0 && s_C >= 0) {
                for(int j=0; j<8; ++j) flux[i-k][j] = get_F_Y(UL, WL, a1L, a2L)[j] + sL*(UL_C[j] - UL[j]);
            } else if(s_C <= 0 && sR >= 0) {
                for(int j=0; j<8; ++j) flux[i-k][j] = get_F_Y(UR, WR, a1R, a2R)[j] + sR*(UR_C[j] - UR[j]);
            } else {
                flux[i-k][0] = UR[2]; 
                flux[i-k][1] = UR[2]*WR.u1;
                flux[i-k][2] = UR[2]*WR.v1 + a1R*WR.P1; 
                flux[i-k][3] = (UR[7] + a1R*WR.P1)*WR.v1;
                flux[i-k][4] = UR[6]; 
                flux[i-k][5] = UR[6]*WR.u2; 
                flux[i-k][6] = UR[6]*WR.v2 + a2R*WR.P2; 
                flux[i-k][7] = (UR[7] + a2R*WR.P2)*WR.v2;
            }
        }
    }
}
