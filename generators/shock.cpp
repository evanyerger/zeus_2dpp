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
; // THIS IS NECESSARY
int main()
{
  std::string sim_dir  = "../../sims/";
  std::string sim_name = "shock";
  std::string path = sim_dir + sim_name + "/";

  struct stat sb;
  if (!(stat(path.c_str(), &sb) == 0) || !S_ISDIR(sb.st_mode))
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
  for (int i=c.nghost; i<c.N1+c.nghost; i++)
    for (int j=c.nghost; j<c.N2+c.nghost; j++)
    {
      if (i < c.N1/2+c.nghost)
      {
        g.d[j*c.full1+i] = 1.0;
        g.e[j*c.full1+i] = 1./(c.gamma - 1);
      }
      else
      {
        g.d[j*c.full1+i] = .125;
        g.e[j*c.full1+i] = .1/(c.gamma - 1);
      }
    }
  for (int i=0; i<c.size; i++) {g.p[i] = (c.gamma-1)*g.e[i];}
  for (int i=0; i<c.size; i++) {g.v1[i] = 0.0;}
  for (int i=0; i<c.size; i++) {g.v2[i] = 0.0;}
  for (int i=0; i<c.size; i++) {g.v3[i] = 0.0;}

  // Apply boundary conditions
  bcs_init();
  update_bcs(&c, &g);

  // Save state file
  int version = 0, ndumps = 0, nsteps = 0;
  double start_time = 0.0, last_dt = 0.0;
  std::string fname = path + "state.txt";
  std::ofstream statefile (fname.c_str());
  statefile << version << std::endl;
  statefile << ndumps  << std::endl;
  statefile << nsteps  << std::endl;
  statefile << start_time << std::endl;
  statefile << last_dt << std::endl;
  statefile << c.N1 << std::endl;
  statefile << c.N2 << std::endl;
  statefile << c.nghost << std::endl;
  statefile.close();

  // Save arrays in grid
  save_variable("p", path, &c, g.p, ndumps);
  save_variable("d", path, &c, g.d, ndumps);
  save_variable("e", path, &c, g.e, ndumps);
#ifdef SELF_GRAVITY
  save_variable("phi", path, &c, g.phi, ndumps);
#endif
  save_variable("v1", path, &c, g.v1, ndumps);
  save_variable("v2", path, &c, g.v2, ndumps);
  save_variable("v3", path, &c, g.v3, ndumps);
#ifdef MHD
  save_variable("b1", path, &c, g.b1, ndumps);
  save_variable("b2", path, &c, g.b2, ndumps);
  save_variable("b3", path, &c, g.b3, ndumps);
  save_variable("emf1", path, &c, g.emf1, ndumps);
  save_variable("emf2", path, &c, g.emf2, ndumps);
  save_variable("emf3", path, &c, g.emf3, ndumps);
#endif

 // free memory
 grid_destruct();
 g.destruct();
}

