#pragma once

#include <iostream>
#include <vector>
#include <algorithm>

#include "../HyperbolicOperator/HyperbolicOperator.h"
#include "../Constants/Physics_Constants.h"
#include "../Constants/Mesh_Constants.h"

class Solver{

private:
    Mesh& mesh;
    const Components phases;
    HyperbolicOperator HyperbolicOp;
public:
    Solver(Mesh& msh);
    void Solve();
private:
    double compute_dt();
};