#include <cmath>
#include "constants.h"
#include "grids.h"

// temporary values
static int fwd, bkwd; //dir=0;

// Pointers to Temporary Arrays
static double *di1=nullptr, *di2=nullptr, *dpar=nullptr,
              *eden=nullptr, *edeni=nullptr,
              *v1i=nullptr, *v2i=nullptr, *v3i=nullptr,
              *s1=nullptr, *s2=nullptr, *s3=nullptr;
static double *D=nullptr, *F=nullptr, *G=nullptr, *H=nullptr, *J=nullptr;

// Private Function Declarations
// Misc
//   specific energy
static void specific_energy(Consts*, Grid*, double* dens);
// Fluxes
//   1-direction
static void interp1(Consts*, Grid*, double*, double*, double);
static void density_flux1(Consts*, Grid*, double*);
static void energy_flux1(Consts*, Grid*, double*);
static void s1_flux1(Consts*, Grid*, double*);
static void s2_flux1(Consts*, Grid*, double*);
static void s3_flux1(Consts*, Grid*, double*);
//   2-direction
static void interp2(Consts*, Grid*, double*, double*, double);
static void density_flux2(Consts*, Grid*, double*);
static void energy_flux2(Consts*, Grid*, double*);
static void s1_flux2(Consts*, Grid*, double*);
static void s2_flux2(Consts*, Grid*, double*);
static void s3_flux2(Consts*, Grid*, double*);
// update variables
//   1-direction
static void update_density1(Consts*, Grid*, double*, double*, double);
static void update_energy1(Consts*, Grid*, double);
static void update_s11(Consts*, Grid*, double*, double);
static void update_s21(Consts*, Grid*, double*, double);
static void update_s31(Consts*, Grid*, double*, double);
//   2-direction
static void update_density2(Consts*, Grid*, double*, double*, double);
static void update_energy2(Consts*, Grid*, double);
static void update_s12(Consts*, Grid*, double*, double);
static void update_s22(Consts*, Grid*, double*, double);
static void update_s32(Consts*, Grid*, double*, double);
// convert momentum to velocity
static void sconvert1(Consts*, Grid*);
static void sconvert2(Consts*, Grid*);
static void sconvert3(Consts*, Grid*);

// Public Functions
//   transport_init(Consts* c, Grid* g)
//     initializes temporary grids for transport step
void transport_init(int size)
{
  di1   = array_init(size);
  di2   = array_init(size);
  dpar  = array_init(size);
  eden  = array_init(size);
  edeni = array_init(size);
  s1    = array_init(size);
  s2    = array_init(size);
  s3    = array_init(size);
  v1i   = array_init(size);
  v2i   = array_init(size);
  v3i   = array_init(size);
  F     = array_init(size);
  G     = array_init(size);
  H     = array_init(size);
  J     = array_init(size);
}

//   transport_step(Consts* c, Grid* g, double dt)
//     main transport step
void transport_step(Consts* c, Grid* g, double dt)
{
  // figure out how to switch directions later

  // update density first
  //   1-direction
  interp1(c, g, g->d, di1, dt); // original density interpolated into dinterp1
  density_flux1(c, g, di1); // density fluxes computed from dinterp1 (interpolated values)
  update_density1(c, g, g->d, dpar, dt);
  //   update specific energy
  specific_energy(c, g, g->d);
  //   2-direction
  interp2(c, g, dpar, di2, dt);
  density_flux2(c, g, di2);

  // update all other variables
  //   1-direction
  interp1(c, g, eden, edeni, dt);
  energy_flux1(c, g, di1);
  update_energy1(c, g, dt);
  interp1(c, g, g->v1, v1i, dt);
  s1_flux1(c, g, di1);
  update_s11(c, g, g->d, dt);
  interp1(c, g, g->v2, v2i, dt);
  s2_flux1(c, g, di1);
  update_s21(c, g, g->d, dt);
  interp1(c, g, g->v3, v3i, dt);
  s3_flux1(c, g, di1);
  update_s31(c, g, g->d, dt);

  update_density2(c, g, dpar, g->d, dt); //WATCH OUT FOR DENSITY FLOOR

  //   2-direction
  specific_energy(c, g, dpar);
  interp2(c, g, eden, edeni, dt);
  energy_flux2(c, g, di2);
  update_energy2(c, g, dt);
  interp2(c, g, g->v1, v1i, dt);
  s1_flux2(c, g, di2);
  update_s12(c, g, g->d, dt);
  interp2(c, g, g->v2, v2i, dt);
  s2_flux2(c, g, di2);
  update_s22(c, g, g->d, dt);
  interp2(c, g, g->v3, v3i, dt);
  s3_flux2(c, g, di2);
  update_s32(c, g, g->d, dt);

  // convert momentum to velocity
  sconvert1(c, g);
  sconvert2(c, g);
  sconvert3(c, g);
}

//   transport_destruct(void)
//     frees memory for temporary arrays
void transport_destruct(void)
{
  array_destruct(di1);
  array_destruct(di2);
  array_destruct(dpar);
  array_destruct(eden);
  array_destruct(edeni);
  array_destruct(s1);
  array_destruct(s2);
  array_destruct(s3);
  array_destruct(v1i);
  array_destruct(v2i);
  array_destruct(v3i);
  array_destruct(F);
  array_destruct(G);
  array_destruct(H);
  array_destruct(J);
}

// Private Functions
//   specific_energy(Consts* c, Grid* g)
//     calculates specific energy
static void specific_energy(Consts* c, Grid* g, double* dens)
{
  for (int i=0; i<c->size; i++)
  {
    eden[i] = g->e[i]/dens[i];
  }
}

//   interp1(Consts* c, double* ary, double dt)
//     interpolates ary using specified method
static void interp1(Consts* c, Grid* g, double* inary, double* outary, double dt)
{
  for (int j=c->nghost-1; j<c->N2+c->nghost+1; j++)
    for (int i=c->nghost-1; i<c->N1+c->nghost+1; i++)
    {
      fwd = j*c->full1 + i; bkwd = fwd - 1;
      outary[fwd] = (g->v1[fwd] > c->v1g) ? inary[bkwd]
#ifdef VLEER
      // valLeero stuff
#endif
      : inary[fwd]
#ifdef VLEER
    // vanLeer stuff
#endif
      ;
    }
}

//   density_flux1(Consts* c, Grid* g, double* dtemp)
//     computes density flux in the 1-direction
static void density_flux1(Consts* c, Grid* g, double* di)
{
  for (int j=c->nghost; j<c->N2+c->nghost+1; j++)
    for (int i=c->nghost; i<c->N1+c->nghost+1; i++)
    {
      fwd = j*c->full1 + i;
      D[fwd] = di[fwd]*(g->v1[fwd]-c->v1g); //leave out dy - divided out later
    }
}

static void update_density1(Consts* c, Grid* g, double* din, double* dout, double dt)
{
  for (int j=c->nghost; j<c->N2+c->nghost; j++)
    for (int i=c->nghost; i<c->N1+c->nghost; i++)
    {
      fwd = j*c->full1 + i + 1; bkwd = fwd - 1;
      dout[bkwd] = din[bkwd]+dt/c->dx*(D[bkwd] - D[fwd]);
    }
}

//   energy_flux1(Consts* c, Grid* g)
//     computes specific energy flux in the 1-direction
static void energy_flux1(Consts* c, Grid* g, double* di)
{
  for (int j=c->nghost; j<c->N2+c->nghost+1; j++)
    for (int i=c->nghost; i<c->N1+c->nghost+1; i++)
    {
      fwd = j*c->full1 + i;
      F[fwd] = edeni[fwd]*di[fwd]*(g->v1[fwd]-c->v1g);
    }
}

//   s1_flux1(Costs* c, Grid* g)
//     calculate fluxes for s1 in 1-direction
static void s1_flux1(Consts* c, Grid* g, double* di)
{
  for (int j=c->nghost-1; j<c->N2+c->nghost; j++)
    for (int i=c->nghost-1; i<c->N1+c->nghost; i++)
    {
      fwd = j*c->full1 + i;
      G[fwd] = .5*v1i[fwd]*(di[fwd]*(g->v1[fwd]-c->v1g) +
          di[fwd+1]*(g->v1[fwd+1]-c->v1g));
    }
}

static void s2_flux1(Consts* c, Grid* g, double* di)
{
  for (int j=c->nghost-1; j<c->N2+c->nghost; j++)
    for (int i=c->nghost-1; i<c->N1+c->nghost; i++)
    {
      fwd = j*c->full1 + i;
      H[fwd] = .5*v2i[fwd]*(di[fwd]*(g->v1[fwd]-c->v1g) +
          di[fwd-c->full1]*(g->v1[fwd-c->full1]-c->v1g));
    }
}

static void s3_flux1(Consts* c, Grid* g, double* di)
{
  for (int j=c->nghost; j<c->N2+c->nghost+1; j++)
    for (int i=c->nghost; i<c->N1+c->nghost+1; i++)
    {
      fwd = j*c->full1 + i;
      J[fwd] = v3i[fwd]*di[fwd]*(g->v1[fwd]-c->v1g);
    }
}

//   interp2(Consts* c, double* ary, double dt)
//     interpolates ary using specified method
static void interp2(Consts* c, Grid* g, double* inary, double* outary, double dt)
{
  for (int j=c->nghost-1; j<c->N2+c->nghost+1; j++)
    for (int i=c->nghost-1; i<c->N1+c->nghost+1; i++)
    {
      fwd = j*c->full1 + i; bkwd = fwd - c->full1;
      outary[fwd] = (g->v2[fwd] > c->v2g) ? inary[bkwd]
#ifdef VLEER
      // valLeero stuff
#endif
      : inary[fwd]
#ifdef VLEER
    // vanLeer stuff
#endif
      ;
    }
}

//   density_flux1(Consts* c, Grid* g, double* dtemp)
//     computes density flux in the 1-direction
static void density_flux2(Consts* c, Grid* g, double* di)
{
  for (int j=c->nghost; j<c->N2+c->nghost+1; j++)
    for (int i=c->nghost; i<c->N1+c->nghost+1; i++)
    {
      fwd = j*c->full1 + i;
      D[fwd] = di[fwd]*(g->v2[fwd]-c->v2g); //leave out dy - divided out later
    }
}

static void update_density2(Consts* c, Grid* g, double* din, double* dout, double dt)
{
  for (int j=c->nghost; j<c->N2+c->nghost; j++)
    for (int i=c->nghost; i<c->N1+c->nghost; i++)
    {
      fwd = (j+1)*c->full1 + i; bkwd = fwd - c->full1;
      dout[bkwd] = din[bkwd]+dt/c->dx*(D[bkwd] - D[fwd]);
    }
}

//   energy_flux1(Consts* c, Grid* g)
//     computes specific energy flux in the 1-direction
static void energy_flux2(Consts* c, Grid* g, double* di)
{
  for (int j=c->nghost; j<c->N2+c->nghost+1; j++)
    for (int i=c->nghost; i<c->N1+c->nghost+1; i++)
    {
      fwd = j*c->full1 + i;
      F[fwd] = edeni[fwd]*di[fwd]*(g->v2[fwd]-c->v2g);
    }
}

//   s1_flux1(Costs* c, Grid* g)
//     calculate fluxes for s1 in 1-direction
static void s1_flux2(Consts* c, Grid* g, double* di)
{
  for (int j=c->nghost-1; j<c->N2+c->nghost; j++)
    for (int i=c->nghost-1; i<c->N1+c->nghost; i++)
    {
      fwd = j*c->full1 + i;
      G[fwd] = .5*v1i[fwd]*(di[fwd]*(g->v1[fwd]-c->v1g) +
          di[fwd+1]*(g->v1[fwd+1]-c->v1g));
    }
}

static void s2_flux2(Consts* c, Grid* g, double* di)
{
  for (int j=c->nghost-1; j<c->N2+c->nghost; j++)
    for (int i=c->nghost-1; i<c->N1+c->nghost; i++)
    {
      fwd = j*c->full1 + i;
      H[fwd] = .5*v2i[fwd]*(di[fwd]*(g->v1[fwd]-c->v1g) +
          di[fwd-c->full1]*(g->v1[fwd-c->full1]-c->v1g));
    }
}

static void s3_flux2(Consts* c, Grid* g, double* di)
{
  for (int j=c->nghost; j<c->N2+c->nghost+1; j++)
    for (int i=c->nghost; i<c->N1+c->nghost+1; i++)
    {
      fwd = j*c->full1 + i;
      J[fwd] = v3i[fwd]*di[fwd]*(g->v1[fwd]-c->v1g);
    }
}

static void update_energy1(Consts* c, Grid* g, double dt)
{
  for (int j=c->nghost; j<c->N2+c->nghost; j++)
    for (int i=c->nghost; i<c->N1+c->nghost; i++)
    {
      fwd = j*c->full1 + i + 1; bkwd = fwd - 1;
      g->e[bkwd] = g->e[bkwd]+dt/c->dx*(F[bkwd] - F[fwd]);
    }
}

static void update_s11(Consts* c, Grid* g, double* di, double dt)
{
  for (int j=c->nghost; j<c->N2+c->nghost; j++)
    for (int i=c->nghost; i<c->N1+c->nghost; i++)
    {
      fwd = j*c->full1 + i; bkwd = fwd - 1;
      s1[fwd] = .5*g->v1[fwd]*(di[fwd]+di[bkwd]) +
        dt/c->dx*(G[bkwd] - G[fwd]);
    }
}

static void update_s21(Consts* c, Grid* g, double* di, double dt)
{
  for (int j=c->nghost; j<c->N2+c->nghost; j++)
    for (int i=c->nghost; i<c->N1+c->nghost; i++)
    {
      fwd = j*c->full1 + i + 1; bkwd = fwd - 1;
      s2[bkwd] = .5*g->v2[bkwd]*(di[bkwd]+di[bkwd-c->full1]) +
        dt/c->dy*(H[bkwd] - H[fwd]);
    }
}

static void update_s31(Consts* c, Grid* g, double* di, double dt)
{
  for (int j=c->nghost; j<c->N2+c->nghost; j++)
    for (int i=c->nghost; i<c->N1+c->nghost; i++)
    {
      fwd = j*c->full1 + i + 1; bkwd = fwd - 1;
      s3[bkwd] = g->v3[bkwd]*di[bkwd] +
        dt/c->dx*(J[bkwd] - J[fwd]);
    }
}

static void update_energy2(Consts* c, Grid* g, double dt)
{
  for (int j=c->nghost; j<c->N2+c->nghost; j++)
    for (int i=c->nghost; i<c->N1+c->nghost; i++)
    {
      fwd = (j+1)*c->full1 + i; bkwd = fwd - c->full1;
      g->e[bkwd] = g->e[bkwd]+dt/c->dx*(F[bkwd] - F[fwd]);
    }
}

static void update_s12(Consts* c, Grid* g, double* di, double dt)
{
  for (int j=c->nghost; j<c->N2+c->nghost; j++)
    for (int i=c->nghost; i<c->N1+c->nghost; i++)
    {
      fwd = (j+1)*c->full1 + i; bkwd = fwd - c->full1;
      s1[fwd] = s1[fwd]*dt/c->dx*(G[bkwd] - G[fwd]);
    }
}

static void update_s22(Consts* c, Grid* g, double* di, double dt)
{
  for (int j=c->nghost; j<c->N2+c->nghost; j++)
    for (int i=c->nghost; i<c->N1+c->nghost; i++)
    {
      fwd = j*c->full1 + i + 1; bkwd = fwd - 1;
      s2[bkwd] = s2[bkwd] + dt/c->dy*(H[bkwd] - H[fwd]);
    }
}

static void update_s32(Consts* c, Grid* g, double* di, double dt)
{
  for (int j=c->nghost; j<c->N2+c->nghost; j++)
    for (int i=c->nghost; i<c->N1+c->nghost; i++)
    {
      fwd = j*c->full1 + i + 1; bkwd = fwd - 1;
      s3[bkwd] = s3[bkwd] + dt/c->dx*(J[bkwd] - J[fwd]);
    }
}

static void sconvert1(Consts* c, Grid* g)
{
  for (int j=c->nghost; j<c->N2+c->nghost; j++)
    for (int i=c->nghost; i<c->N1+c->nghost; i++)
    {
      fwd = j*c->full1 + i; bkwd = fwd - 1;
      g->v1[fwd] = 2*s1[fwd]/(g->d[fwd]+g->d[bkwd]);
    }
}

static void sconvert2(Consts* c, Grid* g)
{
  for (int j=c->nghost; j<c->N2+c->nghost; j++)
    for (int i=c->nghost; i<c->N1+c->nghost; i++)
    {
      fwd = j*c->full1 + i; bkwd = fwd - c->full1;
      g->v2[fwd] = 2*s2[fwd]/(g->d[fwd]+g->d[bkwd]);
    }
}

static void sconvert3(Consts* c, Grid* g)
{
  for (int j=c->nghost; j<c->N2+c->nghost; j++)
    for (int i=c->nghost; i<c->N1+c->nghost; i++)
    {
      fwd = j*c->full1 + i;
      g->v3[fwd] = s3[fwd]/g->d[fwd];
    }
}

