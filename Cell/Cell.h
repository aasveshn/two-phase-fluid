#pragma once

#include <iostream>
#include <string>
#include <vector>
#include "../States/StateW.h"

class Cell{

private:

    unsigned int ID;
    int faces_ID[4];
    
public:

    StateW W;

    Cell(unsigned int ID);
    Cell(unsigned int ID, const int (&faces_id)[4]);

    void setFaces(const  int (&faces_id)[4]);
};