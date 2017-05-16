#include <iostream>
#include <string>
#include <fstream>
#include <ios>
#include <vector>
#include <sys/stat.h>
#include "timestep.h"
#include "constants.h"
#include "grids.h"

static int version;
std::vector<std::string> statevector;
std::vector<double> times;
static int ndumps, nsteps, N1, N2, nghost;
static double last_dt;

void load_variable(std::string, std::string, Consts*, double*, int);
void save_variable(std::string, std::string, Consts*, double*, int);

int check_setup(std::string sim_path)
{
  struct stat sb;
  if (!(stat(sim_path.c_str(), &sb) == 0) || !S_ISDIR(sb.st_mode))
  {
    std::cout << "This simulation has not been set up." << std::endl;
    return 1;
  }
  else {return 0;}
}

void load_sim(std::string path, Consts* c, Grid* g, TimeKeeper* timer, int start_step=-1)
{
  // load state file into vector
  std::ifstream statefile ((path + "state.txt").c_str());
  std::string s;
  while (std::getline(statefile, s))
  {
    statevector.push_back(s);
  }
  statefile.close();
  version = std::stoi(statevector[0]);
  if (version == 0)
  {
    ndumps = std::stoi(statevector[1]);
    nsteps = std::stoi(statevector[2]);
    std::string time_str = statevector[3];
    for (int i=0; i<ndumps+1; i++)
    {
      double t = std::stod(time_str.substr(0, time_str.find(",")));
      times.push_back(t);
    }
    std::cout << times.back() << std::endl;
    last_dt = std::stof(statevector[4]);
    N1 = std::stoi(statevector[5]);
    N2 = std::stoi(statevector[6]);
    nghost = std::stoi(statevector[7]);
  }
  std::cout << ndumps << " " << nsteps << " " << last_dt
    << " " << N1 << " " << N2 << " " << nghost << std::endl;
  // initialize constants
  c->init(N1, N2, nghost);
  // initialize and load grids
  g->init(c->size);
  grid_init(c->size);
  if (start_step < 0) {start_step = ndumps;}
  load_variable("p", path, c, g->p, start_step);
  load_variable("d", path, c, g->d, start_step);
  load_variable("e", path, c, g->e, start_step);
#ifdef SELF_GRAVITY
  load_variable("phi", path, c, g->phi, start_step);
#endif
  load_variable("v1", path, c, g->v1, start_step);
  load_variable("v2", path, c, g->v2, start_step);
  load_variable("v3", path, c, g->v3, start_step);
#ifdef MHD
  load_variable("b1", path, c, g->b1, start_step);
  load_variable("b2", path, c, g->b2, start_step);
  load_variable("b3", path, c, g->b3, start_step);
  load_variable("emf1", path, c, g->emf1, start_step);
  load_variable("emf2", path, c, g->emf2, start_step);
  load_variable("emf3", path, c, g->emf3, start_step);
#endif
  // initialize timer
  timer->init(nsteps, last_dt, times.back());
  if (start_step == 0) {timer->DeltaTZero(c, g);} // set dt to calculated value
}

void save_sim(std::string path)
{
  std::fstream statefile;
  statefile.open((path + "state.txt").c_str(), std::ios::trunc);
  if (!statefile.is_open()) {std::cout << "State is not saved correctly!" << std::endl;}
  statefile << version << std::endl;
  statefile << ndumps  << std::endl;
  statefile << nsteps  << std::endl;
  for (int i=0; i<(int)times.size(); i++) {statefile << times[i] << ",";}
  statefile << std::endl;
  statefile << N1 << std::endl;
  statefile << N2 << std::endl;
  statefile << nghost << std::endl;
  statefile.close();
}

// save all grids into memory and update times
void data_dump(std::string path, Consts* c, Grid* g, double dump_int, TimeKeeper* timer)
{
  if (timer->time > (double)(ndumps + 1)*dump_int)
  {
    ndumps++;
    times.push_back(timer->time);
    last_dt = timer->dt;
    nsteps = timer->nsteps;
    save_variable("p", path, c, g->p, ndumps);
    save_variable("d", path, c, g->d, ndumps);
    save_variable("e", path, c, g->e, ndumps);
#ifdef SELF_GRAVITY
    save_variable("phi", path, c, g->phi, ndumps);
#endif
    save_variable("v1", path, c, g->v1, ndumps);
    save_variable("v2", path, c, g->v2, ndumps);
    save_variable("v3", path, c, g->v3, ndumps);
#ifdef MHD
    save_variable("b1", path, c, g->b1, ndumps);
    save_variable("b2", path, c, g->b2, ndumps);
    save_variable("b3", path, c, g->b3, ndumps);
    save_variable("emf1", path, c, g->emf1, ndumps);
    save_variable("emf2", path, c, g->emf2, ndumps);
    save_variable("emf3", path, c, g->emf3, ndumps);
#endif
  }
  else {;}
}

// load one variable from memory into appropriate grid array
void load_variable(std::string var, std::string path, Consts* c, double* ary, int step=0)
{
  std::string file_number = std::to_string(step);
  size_t fnlength = file_number.length();
  int extra = 4 - fnlength;
  for (int i=0; i<extra; i++) {file_number = "0" + file_number;}
  std::string file_name = var + file_number + ".dat";
  std::ifstream file ((path + var + "/" + file_name).c_str());
  std::string data;
  std::getline(file, data);
  file.close();
  for (int i=0; i<c->size; i++) {ary[i] = std::stod(data.substr(0, data.find(",")));}
}

void save_variable(std::string var, std::string path, Consts* c, double* ary, int step=0)
{
  std::string file_number = std::to_string(step);
  size_t fnlength = file_number.length();
  int extra = 4 - fnlength;
  for (int i=0; i<extra; i++) {file_number = "0" + file_number;}
  std::string file_name = var + file_number + ".dat";
  std::ofstream file ((path + var + "/" + file_name).c_str());
  for (int i=0; i<c->size; i++) {file << ary[i] << ",";}
  file << std::endl;
  file.close();
}
