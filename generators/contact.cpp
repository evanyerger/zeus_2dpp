#include <iostream>
#include <fstream>
#include <ios>
#include <string>
#include <sys/stat.h>
#include "../constants.h"
#include "../grids.h"
#include "../timestep.h"
#include "../bcs.h"
#include "../io.h"

int main()
{
  std::string sim_dir  = "../../sims/";
  std::string sim_name = "sim_name";
  std::string path = sim_dir + sim_name + "/";

  struct stat sb;
  if (!(stat(full_path.c_str(), &sb) == 0) || !S_ISDIR(sb.st_mode))
  {
    std::cout << "The directory required to set up this simulation does not "
      "exist. Please set up with setup.sh" << std::endl;
    return 0;
  }

  // Define simulation constants here
  int N1 = 256, N2 = 256, nghost = 2;
  Consts c;
  c.init(N1, N2, nghost);
  TimeKeeper timer;
  timer.init(0, 0, 0.0);
  Grid g;
  grid_init(c.size);
  g.init(c.size);

  // Set up arrays for p, d, e, phi, v1, v2, v3, b1, b2, b3, etc.


  // Apply boundary conditions
  bcs_init();
  update_bcs(&c, &g);

  // Save state file
  int version = 0, ndumps = 0, nsteps = 0;
  double start_time = 0.0;
  fstream statefile;
  statefile.open((path + "state.txt").c_str(), std::ios::trunc);
  statefile << version << std::endl;
  statefile << ndumps  << std::endl;
  statefile << nsteps  << std::endl;
  statefile << start_time << std::endl;
  statefile << c.N1 << std::endl;
  statefile << c.N2 << std::endl;
  statefile << c.nghost << std::endl;
  statefile.close();

  // Save arrays in grid
  save_variable("p", path, c, g.p, ndumps);
  save_variable("d", path, c, g.d, ndumps);
  save_variable("e", path, c, g.e, ndumps);
#ifdef SELF_GRAVITY
  save_variable("phi", path, c, g.phi, ndumps);
#endif
  save_variable("v1", path, c, g.v1, ndumps);
  save_variable("v2", path, c, g.v2, ndumps);
  save_variable("v3", path, c, g.v3, ndumps);
#ifdef MHD
  save_variable("b1", path, c, g.b1, ndumps);
  save_variable("b2", path, c, g.b2, ndumps);
  save_variable("b3", path, c, g.b3, ndumps);
  save_variable("emf1", path, c, g.emf1, ndumps);
  save_variable("emf2", path, c, g.emf2, ndumps);
  save_variable("emf3", path, c, g.emf3, ndumps);
#endif

 // free memory
 grid_destruct();
 g.destruct();
}

