#ifndef _CONSTANTS_
#define _CONSTANTS_

struct Consts
{
  int N1, N2, nghost, size, full1, full2, start1, end1, start2, end2;
  double dx, dy, dz, dv, gamma, c0, c2, v1g, v2g, v3g;
  void init (int, int, int);
  void SetSpacing(double, double);
};

#endif
