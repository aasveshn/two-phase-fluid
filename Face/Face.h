#pragma once

#include <iostream>

class Face{

public: 

    const unsigned int ID;
    const int Left_ID;
    const int Right_ID;
    const bool isVertical;
    const int type;

    Face(unsigned int id, int L_ID, int R_ID, bool isVert, int Type);

};