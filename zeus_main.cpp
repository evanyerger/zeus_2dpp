#include <iostream>
#include <fstream>
#include "grids.h"
#include "constants.h"
#include "source.h"
#include "timestep.h"

//--------------------------------------------------------------------------//
// Zeus_2d main program
//
// Steps in main:
//   1: Initialize and set constants (set by command line to be added later)
//   2: Initialize grids with appropriate size
//   3: set initial conditions from file
//   4: set function pointers and allociate temporary array memory
//   5: main loop
//   6: final output and memory deallocation
//----------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
  //-1: Initialize and set constants------------------------------------------//
  Consts consts (256, 256, 2);

  //-2: Initialize grids------------------------------------------------------//
  Grid grid (consts.size);

  //-3: Set initial conditions from file and compute initial dt---------------//
  TimeKeeper timer;


  //-4: Allocate temporary memory---------------------------------------------//
  source_init(consts.size);
  grids_init(consts.size);

  //-5: Main Loop-------------------------------------------------------------//
  int maxits = 10000;
  clock_t time = clock();
  while (timer.nits<maxits)
  {
    //---5.1: Source Step-----------------------------------------------------//
#ifdef SELT_GRAVITY
    //poisson solver
#endif
    source_step(&consts, &grid, timer.dt);
    //---5.2: Transport Step--------------------------------------------------//

    //---5.3: Boundary Conditions---------------------------------------------//

    //---5.4: New Timestep----------------------------------------------------//
    timer.DeltaT();
  }
  std::cout << "Main loop took " << clock() - time << " cycles" << std::endl;

  //-6: Final Output and Memory Deallocation----------------------------------//
  source_destruct();
  grids_destruct();
}

/*int main()
{
  // Array Size and Dimensions//
  const int Nx = 256, Ny = 256, ghosts = 2;
  GridSpecs gs (Nx, Ny, ghosts);
  GridSpecs* ptgs = &gs;

  // Define Constants //
  // dx, dy, dt, gamma, C0 (safety factor), C2 (l/dx), vg1, vg2 (grid velocities) //
  int nparams = 8;
  double params[nparams] = {.1, .1, 10, 1.4, .5, 3., 0., 0.};

  // Define Arrays //
  // p, d, e, phi, q1, q2, v1, v2 //
  int nvars = 8;
  double** gridarray = new double* [2*nvars];
  for (int i=0; i<2*nvars; i++) {gridarray[i] = new double[gs.length];}
  int* current_idx = new int[nvars];
  for (int i=0; i<nvars; i++) {current_idx[i] = 0.0;}
  double* temp = new double[gs.length];
  for (int i=0; i<gs.length; i++) {temp[i] = 0.0;}
  Grids grid (ptgs, current_idx, gridarray, params, nvars, temp);
  Grids* ptgrid = &grid;

  // Initial Conditions //
  zeros(ptgs, ptgrid);

  // Main Loop //
  int maxits = 10000;
  for (int it=0; it<maxits; it++)
  {
    // calculate dt //
    grid.deltat();


  }
  return 0;
}*/

