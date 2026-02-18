#pragma once

#include <iostream>
#include <vector>

#include "../States/StateW.h"
#include "../EOS/EOS.h"
#include "../Mesh/Mesh.h"
#include "../Components/Components.h"


class RelaxationOperator{

private:

    Mesh& mesh;
    const Components& phases;

public:

    RelaxationOperator(Mesh& msh, const Components& comp);
    void Relax();
    
private:
    void VelocityRelaxation(Cell& cell);
    void PressureRelaxation(Cell& cell);
    
};