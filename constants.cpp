#include "constants.h"

Consts::Consts(void)
{
  Nx = 128;
  Ny = 128;
  nghost = 2;
  full1 = 2*nghost + Nx;
  full2 = 2*nghost + Ny;
  size = full1*full2;
  start1 = nghost - 1;
  end1 = full1; //?????????????
  start2 = start1;
  end2 = full2; //?????????????
  dx = 1./((double)Nx - 1.);
  dy = 1./((double)Ny - 1.);
  dz = (dx + dy)/2.;
  dv = dx*dy*dz;
  gamma = 5/3;
  c0 = 0.5; //safety factor
  c2 = 3.;
  v1g = 0.;
  v2g = 0.;
  v3g = 0.;
}

Consts::Consts(int nx, int ny, int nGhost)
{
  Nx = nx;
  Ny = ny;
  nghost = nGhost;
  full1 = 2*nghost + Nx;
  full2 = 2*nghost + Ny;
  size = full1*full2;
  start1 = nghost - 1;
  end1 = full1; //?????????????
  start2 = start1;
  end2 = full2; //?????????????
  dx = 1./((double)Nx - 1.);
  dy = 1./((double)Ny - 1.);
  dz = (dx + dy)/2.;
  dv = dx*dy*dz;
  gamma = 5/3;
  c0 = 0.5; //safety factor
  c2 = 3.;
  v1g = 0.;
  v2g = 0.;
  v3g = 0.;
}

void Consts::SetSpacing(double dX, double dY)
{
  dx = dX; dy = dY;
  dz = (dx + dy)/2.;
  dv = dx*dy*dz;
}

