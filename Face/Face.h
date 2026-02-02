#pragma once

#include <iostream>

class Face{

public: 

    const unsigned int Left_ID;
    const unsigned int Right_ID;
    const bool isVertical;

    Face(unsigned int L_ID, unsigned int R_ID, bool isVert);

};