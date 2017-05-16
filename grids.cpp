#include "grids.h"
//#include <iostream>
#include <cmath>

// Pointer to temporary array
double* temp = nullptr;

//_Public functions___________________________________________________________//
//   array_init(int size)
//     Initializes an array
double* array_init(int size)
{
  return new double[size];
}

//   array_destruct(double* ary)
//     Frees an array from memory
void array_destruct(double* ary)
{
  delete[] ary;
}

//   grid_init(int size)
//     allocates memory for temporary array
void grid_init(int size)
{
  temp = array_init(size);
}

//   grid_destruct(void)
//     frees memeory for temporary array
void grid_destruct(void)
{
  array_destruct(temp);
}

//   C_a(Consts* c)
//     Computes adiabatic sound speed
double* C_a(Consts* c, Grid* g)
{
  for (int i=0; i<c->size; i++) {temp[i] = std::sqrt(c->gamma*g->p[i]/g->d[i]);}
  return temp;
}

//   ary_max(Consts* c, double* ary)
//     Finds maximum value in the array
double array_max(Consts* c, double* ary)
{
  double max = ary[0];
  for (int i=1; i<c->size; i++) {max = (ary[i] > max) ? ary[i] : max;}
  return max;
}

//   difference(double* ary, int dir=0)
//     difference array in direction dir
double* difference(Consts* c, double* ary, int dir=0) //difference an array (used in dt calculation)
{
  int fwd, bkwd;
  if (dir == 0)
  {
    for (int j=c->start2; j<c->end2; j++) for (int i=c->start1; i<c->end2; i++)
    {
      fwd = j*c->full2 + i;
      bkwd = j*c->full2 + i - 1;
      temp[fwd] = ary[fwd] - ary[bkwd];
    }
  }
  else
  {
     for (int j=c->start2; j<c->end2; j++) for (int i=c->start1; i<c->end2; i++)
    {
      fwd = j*c->full2 + i;
      bkwd = (j-1)*c->full2 + i;
      temp[fwd] = ary[fwd] - ary[bkwd];
    }
  }
  return temp;
}

// Grids struct
//   Grid::Grid(void)
//     default constructor for Grids struct
Grid::Grid(void)
{
  p = nullptr;
  d = nullptr;
  e = nullptr;
#ifdef SELF_GRAVITY
  phi = nullptr;
#endif
  v1 = nullptr;
  v2 = nullptr;
  v3 = nullptr;
#ifdef MHD
  b1 = nullptr;
  b2 = nullptr;
  b3 = nullptr;
  emf1 = nullptr;
  emf2 = nullptr;
  emf3 = nullptr;
#endif
}


//   Grid::Grid(int size)
//     constructor for Grids struct with array size: size
void Grid::init(int size)
{
  p = array_init(size);
  d = array_init(size);
  e = array_init(size);
#ifdef SELF_GRAVITY
  phi = array_init(size);
#endif
  v1 = array_init(size);
  v2 = array_init(size);
  v3 = array_init(size);
#ifdef MHD
  b1 = array_init(size);
  b2 = array_init(size);
  b3 = array_init(size);
  emf1 = array_init(size);
  emf2 = array_init(size);
  emf3 = array_init(size);
#endif
}

void Grid::destruct(void)
{
  array_destruct(p);
  array_destruct(d);
  array_destruct(e);
#ifdef SELF_GRAVITY
  array_destruct(phi);
#endif
  array_destruct(v1);
  array_destruct(v2);
  array_destruct(v3);
#ifdef MHD
  array_destruct(b1);
  array_destruct(b2);
  array_destruct(b3);
  array_destruct(emf1);
  array_destruct(emf2);
  array_destruct(emf3);
#endif
}

