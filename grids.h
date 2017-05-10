#ifndef _GRIDS_
#define _GRIDS_

#include "constants.h"

// Grids struct declaraiont
struct Grid
{
  double *p, *d, *e, *phi, *v1, *v2, *v3;
#ifdef MHD
  //double *b1, *b2, *b3; OHTERS????
#endif
  Grid ();
  Grid (int);
};

// Public Functions Declarations
double* array_init(int);

void array_destruct(double*);

void grids_init(int);

void grids_destruct(void);

double* C_a(Consts*, Grid*);

double array_max(Consts*, double*);

double* difference(Consts*, double*, int);

#endif

/*#ifndef _GRID_SPECS_
#define _GRID_SPECS_

struct GridSpecs
{
  int Nx, Ny, ghosts, length, fullx, fully, startx, endx, starty, endy;
  GridSpecs(int, int, int);
};

struct Grids
{
  public:
    int nvars;
    double** gridarray;
    int* current_idx;
    double* params;
    GridSpecs* gs;
    double dt, newdt;
    Grids(GridSpecs*, int*, double**, double*, int, double*);
    double* current(int);
    double* next(int);
    double ary_max(double*);
    double* difference(double*, int);
    double avg_c_a();
    deltat();
  private:
    double dt1, dt2, dt3, dt4;
    double* temp;
};

#endif*/
