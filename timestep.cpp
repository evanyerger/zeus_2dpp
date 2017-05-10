#include "constants.h"
#include "grids.h"
#include "timestep.h"
#include <cmath>

// Declare specific variables
double dt1, dt2, dt3, dt4, one, two;

// Public member functions
//   Computes initial dt and configures time-keeping variables starting a
//   simulation
TimeKeeper::Start(Consts* c, Grid* g)
{
  nsteps = 0;
  time = 0.0;
  dt = DeltaTZero(c, g);
  newdt = 0.0;
}

// Restart(int nsteps, double time, double dt)
//   Configures time-keeping variables if restarting a simulation
TimeKeeper::Restart(int a, double b, double c)
{
  nsteps = a;
  time = b;
  dt = c;
  newdt = 0.0;
}

// DeltaT(Consts* c, Grid* g)
//  Computes new dt; updates time and nsteps
void TimeKeeper::DeltaT(Consts* c, Grid* g) //Calculate dt for stability
{
  nsteps++;
  time += dt;
  dt1 = ((c->dx < c->dy) ? c->dx : c->dy)/C_a(c, g);
  dt2 = c->dx / (array_max(g->v1) - c->v1g);
  dt3 = c->dy / (array_max(g->v2) - c->v2g);
  one = .25 * c->dx / c->c2 / array_max(difference(g->v1));
  two = .25 * c->dy / c->c2 / array_max(difference(g->v2));
  dt4 = (one < two) ? one : two;
  newdt = c->c0 / std::sqrt(1./(dt1*dt1) + 1./(dt2*dt2) + 1./(dt3*dt3) + 1./(dt4*dt4));
  dt = (newdt > 1.3*dt) ? 1.3*dt : newdt; // limit new dt to a 30% increase
}

// DeltaTZero(Consts* c, Grid* g)
//  Computes first dt
void TimeKeeper::DeltaT(Consts* c, Grid* g) //Calculate dt for stability
{
  dt1 = ((c->dx < c->dy) ? c->dx : c->dy)/C_a(c, g);
  dt2 = c->dx / (array_max(g->v1) - c->v1g);
  dt3 = c->dy / (array_max(g->v2) - c->v2g);
  one = .25 * c->dx / c->c2 / array_max(difference(g->v1));
  two = .25 * c->dy / c->c2 / array_max(difference(g->v2));
  dt4 = (one < two) ? one : two;
  dt = c->c0 / std::sqrt(1./(dt1*dt1) + 1./(dt2*dt2) + 1./(dt3*dt3) + 1./(dt4*dt4));
}
