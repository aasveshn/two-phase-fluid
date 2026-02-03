#include "Face.h"

Face::Face(unsigned int id, int L_ID, int R_ID, bool isVert, int Type)
            :ID(id), Left_ID(L_ID), Right_ID(R_ID), isVertical(isVert), type(Type){};

