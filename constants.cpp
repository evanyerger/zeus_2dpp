#include "constants.h"

/*Consts::Consts(void)
{
  N1 = 128;
  N2 = 128;
  nghost = 2;
  full1 = 2*nghost + N1;
  full2 = 2*nghost + N2;
  size = full1*full2;
  start1 = nghost;
  end1 = N1 + nghost;
//#ifdef BC_P1
  //start1 = nghost - 1;
  //end1 = N1 + nghost + 1;
//#endif
  start2 = nghost;
  end2 = N2 + nghost;
//#ifdef BC_P2
  //start2 = nghost - 1;
  //end2 = N2 + nghost + 1;
//#endif
  dx = 1./((double)N1 - 1.);
  dy = 1./((double)N2 - 1.);
  dz = (dx + dy)/2.;
  dv = dx*dy*dz;
  gamma = 5/3;
  c0 = 0.5; //safety factor
  c2 = 3.;
  v1g = 0.;
  v2g = 0.;
  v3g = 0.;
}*/

void Consts::init(int nx, int ny, int nGhost)
{
  N1 = nx;
  N2 = ny;
  nghost = nGhost;
  full1 = 2*nghost + N1;
  full2 = 2*nghost + N2;
  size = full1*full2;
  start1 = nghost;
  end1 = N1 + nghost;
//#ifdef BC_P1
  //start1 = nghost - 1;
  //end1 = N1 + nghost + 1;
//#endif
  start2 = nghost;
  end2 = N2 + nghost;
//#ifdef BC_P2
  //start2 = nghost - 1;
  //end2 = N2 + nghost + 1;
//#endif
  dx = 1./((double)N1 - 1.);
  dy = 1./((double)N2 - 1.);
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

