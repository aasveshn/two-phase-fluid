#pragma once

#include <iostream>
#include <vector>
#include <algorithm>

#include "../States/StateW.h"
#include "../States/StateU.h"
#include "../Face/Face.h"
#include "../Cell/Cell.h"
#include "../Mesh/Mesh.h"
#include "../Components/Components.h"
#include "../EOS/EOS.h"

class RiemanSolver{

private:

    std::vector<StateU> Flux_X;
    std::vector<StateU> Flux_Y;

public:

    RiemanSolver(unsigned int N, unsigned int k);
    void HLLC(const Mesh& mesh, const std::vector<StateW>& WL_in, const std::vector<StateW>& WR_in,
              unsigned int H_idx, bool is_X_dir, const Components& phases);
    const std::vector<StateU>& getXFlux() const {return Flux_X;};
    const std::vector<StateU>& getYFlux() const {return Flux_Y;};

};