#pragma once

#include <iostream>
#include <vector>
#include <algorithm>

#include "../HyperbolicOperator/HyperbolicOperator.h"
#include "../Constants/Physics_Constants.h"
#include "../Constants/Mesh_Constants.h"
#include "../RelaxationOperator/RelaxationOperator.h"

class Solver{

private:
    Mesh& mesh;
    const Components phases;
    HyperbolicOperator HyperbolicOp;
    RelaxationOperator RelaxationOp;

public:
    Solver(Mesh& msh);
    const Components& getPhases();
    void Solve();
private:
    double compute_dt();
};