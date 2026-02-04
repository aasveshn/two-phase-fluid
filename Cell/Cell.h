#pragma once

#include <iostream>
#include <string>
#include <vector>
#include "../States/StateW.h"

class Cell{

public:

    const unsigned int ID, i, j;
    //0 - левая
    //1 - правая
    //2 - верхняя
    //3 - нижняя
    const unsigned int faces_ID[4];
    

    StateW W;

    
    Cell(unsigned int id, unsigned int i_index, unsigned int j_index, const unsigned int (&faces_id)[4]);
    Cell(const Cell& other);
    
};