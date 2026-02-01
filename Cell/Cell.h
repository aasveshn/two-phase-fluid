#pragma once

#include <iostream>
#include <string>
#include <vector>
#include "../States/StateW.h"

class Cell{

private:

    unsigned int ID;
    unsigned int faces_ID[4];
    
public:

    StateW W;

    Cell(unsigned int ID);
    Cell(unsigned int ID, const unsigned int (&faces_id)[4]);

    void setFaces(const unsigned int (&faces_id)[4]);
};