#ifndef _TIMESTEP_
#define _TIMESTEP_

#include "constants.h"
#include "grids.h"

void timeKeeper_init(int);

void timeKeeper_destruct(int);

struct TimeKeeper
{
  int nsteps;
  double dt, newdt, time;
  void Start(Consts*, Grid*);
  void Restart(int, double, double);
  void DeltaT(Consts* , Grid* );
};

#endif
