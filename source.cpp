#include <iostream>
#include "constants.h"
#include "grids.h"
#include "source.h"

// Temporary Values
static int fwd, bkwd, bkwd2;
double val, divv;

// Pointers to temporary arrays
static double *v1na=nullptr, *v2na=nullptr, *v1nb=nullptr, *v2nb=nullptr;
static double *q1=nullptr, *q2=nullptr, *enb=nullptr;

// Private function declarations
static void pressure(Consts*, Grid*, double);
static void subOneX(Consts*, Grid*, double);
static void subOneY(Consts*, Grid*, double);
static void subTwoQ(Consts* c, Grid* g, double);
static void subTwoVVE(Consts* c, Grid* g, double);
static void subThree(Consts* c, Grid* g, double);
static void replaceV(Consts* c, Grid* g);

//---PUBLIC FUNCTIONS---------------------------------------------------------//
//   source_init(int size):
//     allocates memory for temporary arrays and sets their pointers

void source_init(int size)
{
  v1na = array_init(size);
  v2na = array_init(size);
  v1nb = array_init(size);
  v2nb = array_init(size);
  q1   = array_init(size);
  q2   = array_init(size);
  enb  = array_init(size);
}

//   source_step(Consts* c, Grid* g)
//     performs source step

void source_step(Consts* c, Grid* g, double dt)
{
  //   calculate pressure from internal energy
  pressure(c, g, c->gamma-1);

  //   First substep
  subOneX(c, g, dt);
  subOneY(c, g, dt);

  //   Second substep
  subTwoQ(c, g, dt);
  subTwoVVE(c, g, dt);
  //   Third substep
  subThree(c, g, dt);
  //   Replace v1, v2
  replaceV(c, g);
}

//   source_destruct()
//     frees memory allocated for temporary arrays
void source_destruct(void)
{
  array_destruct(v1na);
  array_destruct(v2na);
  array_destruct(v1nb);
  array_destruct(v2nb);
  array_destruct(q1);
  array_destruct(q2);
  array_destruct(enb);
}

//---PRIVATE FUNCTIONS--------------------------------------------------------//
//   pressure(Consts* c, Grid* g, double gamma)
//     computes pressure from internal energy using ideal equation of state
static void pressure(Consts* c, Grid* g, double gamone)
{
  for (int i=0; i<c->size; i++) g->p[i] = gamone*g->e[i];
}

//   subOneX(Consts* c, Grid* g, double dt)
//     source velocity change from pressure gradient and self-gravity in the
//     1-direction
static void subOneX(Consts* c, Grid* g, double dt)
{
  for (int j=c->start2; j<c->end2; j++) for (int i=c->start1; i<c->end1; i++)
  {
    fwd = j*c->full1 + i; bkwd = j*c->full1 + i - 1;

    v1na[fwd] = g->v1[fwd] - 2.*dt/c->dx*((g->p[fwd] - g->p[bkwd])/
      (g->d[fwd] + g->d[bkwd])
#ifdef SELF_GRAVITY
      + g->phi[fwd] - g->phi[bkwd]
#endif
      );
  }
}

//   subOneY(Consts* c, Grid* g, double dt)
//     source velocity change from pressure gradient and self-gravity in the
//     2-direction
static void subOneY(Consts* c, Grid* g, double dt)
{
  for (int j=c->start2; j<c->end2; j++) for (int i=c->start1; i<c->end1; i++)
  {
    fwd = j*c->full1 + i; bkwd = (j-1)*c->full1 + i;
    v2na[fwd] = g->v2[fwd] - 2.*dt/c->dy*((g->p[fwd] - g->p[bkwd])/
      (g->d[fwd] + g->d[bkwd])
#ifdef SELF_GRAVITY
      + g->phi[fwd] - g->phi[bkwd]
#endif
      );
  }
}

//   subTwoQ(Consts* c, Grid* g, double dt)
//     compute viscous pressure in 1- and 2-directions
static void subTwoQ(Consts* c, Grid* g, double dt)
{
  for (int j=c->start2; j<c->end2; j++) for (int i=c->start1; i<c->end1; i++)
  {
    bkwd = j*c->full1 + i;
    val = g->v1[j*c->full1 + i + 1] - g->v1[bkwd];
    q1[bkwd] = (val < 0) ? c->c2*g->d[bkwd]*val*val : 0.0;
    val = g->v2[(j+1)*c->full1 + i] - g->v2[bkwd];
    q2[bkwd] = (val < 0) ? c->c2*g->d[bkwd]*val*val : 0.0;
  }
}

//   subTwoVVE(Consts* c, Grid* g, double dt)
//     update velocity from viscous pressure and add viscous heating
static void subTwoVVE(Consts* c, Grid* g, double dt)
{
  for (int j=c->start2; j<c->end2; j++) for (int i=c->start1; i<c->end1; i++)
  {
    fwd = j*c->full1 + i; bkwd = j*c->full1 + i - 1; bkwd2 = (j-1)*c->full1 + i;
    v1nb[fwd] = v1na[fwd] + 2*dt/c->dx*(q1[bkwd] - q1[fwd])/(g->d[fwd] + g->d[bkwd]);
    v2nb[fwd] = v2na[fwd] + 2*dt/c->dy*(q1[bkwd2] - q1[fwd])/(g->d[fwd] + g->d[bkwd2]);
    enb[fwd] = g->e[fwd] - dt*(q1[fwd]*(g->v1[bkwd+2] - g->v1[fwd])/c->dx +
        q2[fwd]*(g->v2[bkwd2+2] - g->v2[fwd])/c->dy);
  }
}

//   subThree(Consts* c, Grid* g, double dt)
//     add compressional heating
static void subThree(Consts* c, Grid* g, double dt)
{
  val = c->gamma - 1;
  for (int j=c->start2; j<c->end2; j++) for (int i=c->start1; i<c->end1; i++)
  {
    fwd = j*c->full1 + i;
    divv = (g->v1[fwd+1] - g->v1[fwd])/c->dx + (g->v2[(j+1)*c->full1+i] - g->v2[fwd])/c->dy;
    g->e[fwd] = (1. - .5*dt*val*divv)/(1 + .5*dt*val*divv)*enb[fwd];
  }
}

//   replaceV(Consts* c, Grid* g)
//     replaces V1 and V2 with V1nb and V2nb, respectively
static void replaceV(Consts* c, Grid* g)
{
  for (int i=0; i<c->size; i++)
  {
    g->v1[i] = v1nb[i]; g->v2[i] = v2nb[i];
  }
}


