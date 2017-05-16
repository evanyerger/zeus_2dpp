#ifndef _TIMESTEP_
#define _TIMESTEP_

#include "constants.h"
#include "grids.h"

struct TimeKeeper
{
  int nsteps;
  double dt, newdt, time;
  void init(int, double, double);
  void DeltaT(Consts*, Grid*);
  void DeltaTZero(Consts*, Grid*);
};

#endif
