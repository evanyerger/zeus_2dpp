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

//   grids_init(int size)
//     Initializes the temporary array(s) needed for public functions
void grids_init(int size)
{
  temp = array_init(size);
}

//   grids_destruct(void)
//     Frees memory ossociated with temporary array(s) used in public functions
void grids_destruct(void)
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
  p = array_init(128*128);
  d = array_init(128*128);
  e = array_init(128*128);
#ifdef SELF_GRAVITY
  phi = array_init(128*128);
#endif
  v1 = array_init(128*128);
  v2 = array_init(128*128);
#ifdef MHD
  b1 = array_init(128*128);
  b2 = array_init(128*128);
  b3 = array_init(128*128);
#endif
}


//   Grid::Grid(int size)
//     constructor for Grids struct with array size: size
Grid::Grid(int size)
{
  p = array_init(size);
  d = array_init(size);
  e = array_init(size);
#ifdef SELF_GRAVITY
  phi = array_init(size);
#endif
  v1 = array_init(size);
  v2 = array_init(size);
#ifdef MHD
  b1 = array_init(size);
  b2 = array_init(size);
  b3 = array_init(size);
#endif
}


/* Container for information on the grids */

/*GridSpecs::GridSpecs(int a, int b, int c)
{
  Nx = a;
  Ny = b;
  ghosts = c;
  fullx = Nx + 2*ghosts;
  fully = Ny + 2*ghosts;
  length = fullx * fully;
  startx = ghosts;
  starty = ghosts;
  endx = Nx - ghosts;
  endy = Ny - ghosts;
}

Grids::Grids(GridSpecs* ptgs, int* a, double** grida, double* b, int c, double* t)
{
  gs = ptgs;
  current_idx = a;
  params = b;
  nvars = c;
  gridarray = grida;
  temp = t;
  dt1 = 0;
  dt2 = 0;
  dt3 = 0;
  dt4 = 0;
  newdt = 0;
  dt = 1000;
}

double* Grids::current(int i) //returns current array of specified variable
{
  return gridarray[2*i + current_idx[i]];
}

double* Grids::next(int i) //returns free array of specified variable
{
  current_idx[i] = (current_idx[i] + 1) % 2;
  return gridarray[2*i + current_idx[i]];
}

double Grids::ary_max(double* ary) //find maximum of an array
{
  double max = ary[0];
  for (int i=1; i<gs->length; i++) {if (ary[i] > max) {max = ary[i];}}
  return max;
}

double* Grids::difference(double* ary, int dir=0) //difference an array (used in dt calculation)
{
  int fwd, bkwd;
  if (dir == 0)
  {
    for (int i=gs->startx; i<gs->endx; i++) for (int j=starty; j<gs->endy; j++)
    {
      fwd = i*gs->fullx + j;
      bkwd = (i - 1)*gs->fullx + j;
      temp[fwd] = ary[fwd] - ary[bkwd];
    }
  }
  else
  {
    for (int i=gs->startx; i<gs->endx; i++) for (int j=starty; j<gs->endy; j++)
    {
      fwd = i*gs->fullx + j;
      bkwd = i*gs->fullx + j - 1;
      temp[fwd] = ary[fwd] - ary[bkwd];
    }
  }
  return &temp;
}

double Grids::avg_c_a() //Average adiabatic sound speed
{
  double sum = 0.0;
  for (int i=0; i<gs->length; i++) {sum += current(0)[i]/current(1)[i];}
  return params[3] * sum;
}

void Grids::deltat() //Calculate dt for stability
{
  dt1 = ((params[0] < params[1]) ? params[0] : params[1])/avg_c_a();
  dt2 = params[0] / (ary_max(current(6)) - params[6]);
  dt3 = params[1] / (ary_max(current(7)) - params[7]);
  double one = .25 * params[0] / params[5] / ary_max(difference(current(6)));
  double two = .25 * params[1] / params[5] / ary_max(difference(current(7)));
  dt4 = (one < two) ? one : two;
  newdt = params[4] / std::sqrt(1./(dt1*dt1) + 1./(dt2*dt2) + 1./(dt3*dt3) + 1./(dt4*dt4));
  dt = (newdt > 1.3*dt) ? 1.3*dt : newdt;
}*/
