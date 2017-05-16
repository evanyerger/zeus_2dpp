#include "constants.h"
#include "grids.h"
#include "timestep.h"
#include <cmath>

// Declare specific variables
double dt1, dt2, dt3, dt4, one, two;

// Public member functions
// init(int nsteps, double time, double dt)
//   Configures time-keeping variables from given values
void TimeKeeper::init(int a, double b, double c)
{
  nsteps = a;
  dt = b;
  time = c;
}

// DeltaT(Consts* c, Grid* g)
//  Computes new dt; updates time and nsteps
void TimeKeeper::DeltaT(Consts* c, Grid* g) //Calculate dt for stability
{
  nsteps++;
  time += dt;
  dt1 = ((c->dx < c->dy) ? c->dx : c->dy)/array_max(c, C_a(c, g));
  dt2 = c->dx / (array_max(c, g->v1) - c->v1g);
  dt3 = c->dy / (array_max(c, g->v2) - c->v2g);
  one = .25 * c->dx / c->c2 / array_max(c, difference(c, g->v1, 1));
  two = .25 * c->dy / c->c2 / array_max(c, difference(c, g->v2, 2));
  dt4 = (one < two) ? one : two;
  newdt = c->c0 / std::sqrt(1./(dt1*dt1) + 1./(dt2*dt2) + 1./(dt3*dt3) + 1./(dt4*dt4));
  dt = (newdt > 1.3*dt) ? 1.3*dt : newdt; // limit new dt to a 30% increase
}

// DeltaTZero(Consts* c, Grid* g)
//  Computes first dt
void TimeKeeper::DeltaTZero(Consts* c, Grid* g) //Calculate dt for stability
{
  dt1 = ((c->dx < c->dy) ? c->dx : c->dy)/array_max(c, C_a(c, g));
  dt2 = c->dx / (array_max(c, g->v1) - c->v1g);
  dt3 = c->dy / (array_max(c, g->v2) - c->v2g);
  one = .25 * c->dx / c->c2 / array_max(c, difference(c, g->v1, 1));
  two = .25 * c->dy / c->c2 / array_max(c, difference(c, g->v2, 2));
  dt4 = (one < two) ? one : two;
  dt = c->c0 / std::sqrt(1./(dt1*dt1) + 1./(dt2*dt2) + 1./(dt3*dt3) + 1./(dt4*dt4));
}
