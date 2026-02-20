#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>

#include "../Cell/Cell.h"
#include "../Face/Face.h"
#include "../EOS/EOS.h"
#include "../Components/Components.h"

class Mesh{

public:

    //const double hx, hy;
    std::vector<Cell> Cells;
    std::vector<Face> Faces;
    unsigned int VertToHoriz;

    Mesh(std::string filename);
    void SetInitialCondidions(const Components& phases);
};