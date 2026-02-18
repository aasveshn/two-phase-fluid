#pragma once

#include <iostream>
#include <vector>

#include "../Mesh/Mesh.h"
#include "../Components/Components.h"
#include "../MUSCL/MUSCL.h"
#include "../RiemanSolver/RiemanSolver.h"
#include "../States/StateU.h"
#include "../States/StateW.h"

class HyperbolicOperator{

private:

    const Components& phases;
    MUSCL muscl;
    RiemanSolver rieamn;
    Mesh& mesh;

public:

    HyperbolicOperator(Mesh& msh, const Components& comp);
    void HyperbolicStepX(double dt);
    void HyperbolicStepY(double dt);

private:
    void GodunovStep(bool is_X_dir, double dt);

};