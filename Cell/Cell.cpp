#include "Cell.h"



Cell::Cell(unsigned int id,
           unsigned int i_index,
           unsigned int j_index,
           const unsigned int (&faces_id)[4])
           : ID(id),
           i(i_index),
           j(j_index),
           faces_ID{ faces_id[0], faces_id[1], faces_id[2], faces_id[3] }{}


Cell::Cell(const Cell& other) 
    : ID(other.ID), 
      i(other.i), 
      j(other.j), 
      faces_ID{other.faces_ID[0], other.faces_ID[1], other.faces_ID[2], other.faces_ID[3]},
      W(other.W){}

bool Cell::is_Interface()
{
  return W.a1 > 1e-4 && W.a1 < 1-1e-4;
}