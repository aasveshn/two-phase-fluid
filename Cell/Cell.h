#pragma once

#include <iostream>
#include <string>
#include <vector>
#include "../States/StateW.h"

class Cell{

public:

    const unsigned int ID, i, j;
    const unsigned int faces_ID[4];
    

    StateW W;

    //Cell(unsigned int ID);
    Cell(unsigned int id, unsigned int i_index, unsigned int j_index, const unsigned int (&faces_id)[4]);

    //void setFaces(const  int (&faces_id)[4]);
};