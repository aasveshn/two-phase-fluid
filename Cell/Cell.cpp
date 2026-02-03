#include "Cell.h"

Cell::Cell(unsigned int id){ ID = id; }
Cell::Cell(unsigned int id, const int (&faces_id)[4])
{
    ID = id; 
    faces_ID[0] = faces_id[0]; 
    faces_ID[1] = faces_id[1];
    faces_ID[2] = faces_id[2];
    faces_ID[3] = faces_id[3];
}

void Cell::setFaces(const int (&faces_id)[4])
{
    faces_ID[0] = faces_id[0]; 
    faces_ID[1] = faces_id[1];
    faces_ID[2] = faces_id[2];
    faces_ID[3] = faces_id[3];
}

