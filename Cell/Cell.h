#pragma once

#include <iostream>
#include <string>
#include <vector>
#include "../State/State.h"

class Cell{

private:

    unsigned int ID;
    unsigned int faces_ID[4];
    
public:

    State W;

    Cell(unsigned int ID);
    Cell(unsigned int ID, const unsigned int (&faces_id)[4]);

    void setFaces(const unsigned int (&faces_id)[4]);
    std::vector<double> getU();
};