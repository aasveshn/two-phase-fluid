#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

#include "../Cell/Cell.h"
#include "../Face/Face.h"

class Mesh{

public:

    const double hx, hy;
    

    std::vector<Cell> Cells;
    std::vector<Face> Faces;

    
};