#pragma once

#include <iostream>
#include <vector>
#include <algorithm>

#include "../EOS/EOS.h"
#include "../Mesh/Mesh.h"
#include "../States/StateW.h"
#include "../Constants/Mesh_Constants.h"
#include "../Components/Components.h"

class MUSCL{

private:

    const double B;
    std::vector<StateW> W_L;
    std::vector<StateW> W_R;
    //std::vector<StateW> W_ex_L;
    //std::vector<StateW> W_ex_R;

public:
    
    MUSCL(int N);
    void MUSCL_Operator(const Mesh& mesh, bool is_X_dir, double dt, const Components& phases);
    const std::vector<StateW>& getW_L() {return W_L;}
    const std::vector<StateW>& getW_R() {return W_R;}

private:

    StateW slopLimiter(const StateW& dL, const StateW& dR);
    StateW A_prod_W(const StateW& W, const StateW& WL, const StateW& WR, const Components& phases);
};

