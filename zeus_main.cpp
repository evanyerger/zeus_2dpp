#include <iostream>
#include <fstream>
#include <string>
#include "io.h"
#include "grids.h"
#include "constants.h"
#include "source.h"
#include "timestep.h"
#include "bcs.h"
#include "transport.h"

static void print_variable(Consts* c, double* ary)
{
 for (int j=c->full2-1; j>=0; j--)
 {
   for (int i=0; i<c->full1; i++)
   {
     std::cout << ary[j*c->full1 + i] << " ";
   }
 std::cout << std::endl;
  }
}
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
  //-1: Initialize, set constants, and load initial conditions----------------//
  std::string sim_dir = "../sims/";
  if (argc < 2)
  {
    std::cout << "Please provide a simulation name." << std::endl;
    return 0;
  }
  std::string sim_name = std::string(argv[1]);
  std::string path = sim_dir + sim_name + "/";
  int check = check_setup(path);
  if (check == 1) {return 0;}
  Consts c;
  Grid g;
  TimeKeeper timer;
  std::cout << "start load" << std::endl;
  load_sim(path, &c, &g, &timer, 0);
  std::cout << "dt0: " << timer.dt << std::endl;
  //-4: Set function pointers and allocate temporary memory-------------------//
  source_init(c.size);
  transport_init(c.size);
  bcs_init();

  //-5: Main Loop-------------------------------------------------------------//
  std::cout << "Starting Main Loop" << std::endl;
  double max_time = 0.101, dump_int = .01;
  clock_t time = clock();
  while (timer.time<max_time)
  {
    //---5.1: Source Step-----------------------------------------------------//
#ifdef SELF_GRAVITY
    //poisson solver
#endif
    source_step(&c, &g, timer.dt);
    update_bcs(&c, &g);
    //---5.2: Transport Step--------------------------------------------------//
    transport_step(&c, &g, timer.dt);
    //---5.3: Boundary Conditions---------------------------------------------//
    update_bcs(&c, &g);

    //---5.4: Data IO---------------------------------------------------------//
    data_dump(path, &c, &g, dump_int, &timer);

    //---5.4: New Timestep----------------------------------------------------//
    timer.DeltaT(&c, &g);
  }
  std::cout << "Main loop took " << clock() - time << " cycles" << std::endl;

  //-6: Final Output and Memory Deallocation----------------------------------//
  save_sim(path, timer.dt);
  g.destruct();
  grid_destruct();
  source_destruct();
  transport_destruct();
}

