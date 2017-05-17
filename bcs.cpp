#include <iostream>
#include "constants.h"
#include "grids.h"

int idx, ghstidx;

// function pointers to boundary conditions updators
// Left, Right, Top, and Bottom
void (*bcL)(Consts*, Grid*) = nullptr;
void (*bcR)(Consts*, Grid*) = nullptr;
void (*bcT)(Consts*, Grid*) = nullptr;
void (*bcB)(Consts*, Grid*) = nullptr;
void (*bL)(Consts*, double*) = nullptr;
void (*bR)(Consts*, double*) = nullptr;
void (*bT)(Consts*, double*) = nullptr;
void (*bB)(Consts*, double*) = nullptr;


// forward declaration of private functions
#ifdef BC_P1
static void periodicL(Consts*, Grid*);
static void periodicR(Consts*, Grid*);
static void pL(Consts*, double*);
static void pR(Consts*, double*);
#endif
#ifdef BC_P1
static void periodicT(Consts*, Grid*);
static void periodicB(Consts*, Grid*);
static void pT(Consts*, double*);
static void pB(Consts*, double*);
#endif

//-Public Functions-----------------------------------------------------------//
// bcs_init(void)
//   sets appropriate pointers to private bc functions
void bcs_init(void)
{
#ifdef BC_P1
  bcL = &periodicL;
  bcR = &periodicR;
  bL = &pL;
  bR = &pR;
#endif
#ifdef BC_P2
  bcT = &periodicT;
  bcB = &periodicB;
  bT = &pT;
  bB = &pB;
#endif
}

// update_bcs(Consts* c, Grid* g)
//   updates boundary conditions on appropriate variables
void update_bcs(Consts* c, Grid* g)
{
  (*bcL) (c, g);
  (*bcR) (c, g);
  (*bcT) (c, g);
  (*bcB) (c, g);
  // update corners
  for (int i=0; i<c->nghost; i++) for (int j=0; j<c->nghost; j++)
  {
    g->p[j*c->full1 + i] = g->p[(j+c->N2)*c->full1 + i + c->N1];
    g->p[j*c->full1 + i + c->N1 + c->nghost] = g->p[(j+c->N2)*c->full1 + i + c->nghost];
    g->p[(j+c->N2+c->nghost)*c->full1 + i] = g->p[(j+c->nghost)*c->full1 + i + c->N1];
    g->p[(j+c->N2+c->nghost)*c->full1 +i+c->N1+c->nghost] = g->p[(j+c->nghost)*c->full1 +i+c->nghost];
    g->d[j*c->full1 + i] = g->d[(j+c->N2)*c->full1 + i + c->N1];
    g->d[j*c->full1 + i + c->N1 + c->nghost] = g->d[(j+c->N2)*c->full1 + i + c->nghost];
    g->d[(j+c->N2+c->nghost)*c->full1 + i] = g->d[(j+c->nghost)*c->full1 + i + c->N1];
    g->d[(j+c->N2+c->nghost)*c->full1 +i+c->N1+c->nghost] = g->d[(j+c->nghost)*c->full1 +i+c->nghost];
    g->e[j*c->full1 + i] = g->e[(j+c->N2)*c->full1 + i + c->N1];
    g->e[j*c->full1 + i + c->N1 + c->nghost] = g->e[(j+c->N2)*c->full1 + i + c->nghost];
    g->e[(j+c->N2+c->nghost)*c->full1 + i] = g->e[(j+c->nghost)*c->full1 + i + c->N1];
    g->e[(j+c->N2+c->nghost)*c->full1 +i+c->N1+c->nghost] = g->e[(j+c->nghost)*c->full1 +i+c->nghost];
#ifdef SELF_GRAVITY
    g->phi[j*c->full1 + i] = g->phi[(j+c->N2)*c->full1 + i + c->N1];
    g->phi[j*c->full1 + i + c->N1 + c->nghost] = g->phi[(j+c->N2)*c->full1 + i + c->nghost];
    g->phi[(j+c->N2+c->nghost)*c->full1 + i] = g->phi[(j+c->nghost)*c->full1 + i + c->N1];
    g->phi[(j+c->N2+c->nghost)*c->full1 +i+c->N1+c->nghost] = g->phi[(j+c->nghost)*c->full1 +i+c->nghost];
#endif
    g->v1[j*c->full1 + i] = g->v1[(j+c->N2)*c->full1 + i + c->N1];
    g->v1[j*c->full1 + i + c->N1 + c->nghost] = g->v1[(j+c->N2)*c->full1 + i + c->nghost];
    g->v1[(j+c->N2+c->nghost)*c->full1 + i] = g->v1[(j+c->nghost)*c->full1 + i + c->N1];
    g->v1[(j+c->N2+c->nghost)*c->full1 +i+c->N1+c->nghost] = g->v1[(j+c->nghost)*c->full1 +i+c->nghost];
    g->v2[j*c->full1 + i] = g->v2[(j+c->N2)*c->full1 + i + c->N1];
    g->v2[j*c->full1 + i + c->N1 + c->nghost] = g->v2[(j+c->N2)*c->full1 + i + c->nghost];
    g->v2[(j+c->N2+c->nghost)*c->full1 + i] = g->v2[(j+c->nghost)*c->full1 + i + c->N1];
    g->v2[(j+c->N2+c->nghost)*c->full1 +i+c->N1+c->nghost] = g->v2[(j+c->nghost)*c->full1 +i+c->nghost];
    g->v3[j*c->full1 + i] = g->v3[(j+c->N2)*c->full1 + i + c->N1];
    g->v3[j*c->full1 + i + c->N1 + c->nghost] = g->v3[(j+c->N2)*c->full1 + i + c->nghost];
    g->v3[(j+c->N2+c->nghost)*c->full1 + i] = g->v3[(j+c->nghost)*c->full1 + i + c->N1];
    g->v3[(j+c->N2+c->nghost)*c->full1 +i+c->N1+c->nghost] = g->v3[(j+c->nghost)*c->full1 +i+c->nghost];
  }
}

void single_bc(Consts* c, double* ary)
{
  (*bL) (c, ary);
  (*bR) (c, ary);
  (*bT) (c, ary);
  (*bB) (c, ary);
  // update corners
  for (int i=0; i<c->nghost; i++) for (int j=0; j<c->nghost; j++)
  {
    ary[j*c->full1 + i] = ary[(j+c->N2)*c->full1 + i + c->N1];
    ary[j*c->full1 + i + c->N1 + c->nghost] = ary[(j+c->N2)*c->full1 + i + c->nghost];
    ary[(j+c->N2+c->nghost)*c->full1 + i] = ary[(j+c->nghost)*c->full1 + i + c->N1];
    ary[(j+c->N2+c->nghost)*c->full1 +i+c->N1+c->nghost] = ary[(j+c->nghost)*c->full1 +i+c->nghost];
  }
}

//-Private Functions----------------------------------------------------------//
#ifdef BC_P1
static void periodicL(Consts* c, Grid* g)
{
  // start1 and end1 are shifted by 1 for periodic conditions
  // -- work from primary variables
  for (int j=c->nghost; j<c->N2+c->nghost; j++) for (int i=0; i<c->nghost; i++)
  {
    ghstidx = j*c->full1 + i; idx = ghstidx + c->N1;
    g->p[ghstidx] = g->p[idx];
    g->d[ghstidx] = g->d[idx];
    g->e[ghstidx] = g->e[idx];
#ifdef SEFL_GRAVITY
    g->phi[ghstidx] = g->phi[idx];
#endif
    g->v1[ghstidx] = g->v1[idx];
    g->v2[ghstidx] = g->v2[idx];
    g->v3[ghstidx] = g->v3[idx];
#ifdef MHD
    g->emf1[ghstidx] = g->emf1[idx];
    g->emf2[ghstidx] = g->emf2[idx];
    g->emf3[ghstidx] = g->emf3[idx];
#endif
  }
  // normal velocity
  /*for (int j=c->nghost; j<c->nghost+c->N2; j++)
  {
    ghstidx = j*c->full2;
    g->v1[ghstidx] = g->v1[ghstidx + c->N1];
  }*/
}

static void periodicR(Consts* c, Grid* g)
{
  // start and end variables shifted by 1 for periodic conditions
  // -- work from primary variables
  for (int j=c->nghost; j<c->N2+c->nghost; j++) for (int i=c->N1+c->nghost; i<c->full1; i++)
  {
    ghstidx = j*c->full1 + i; idx = ghstidx - c->N1;
    g->p[ghstidx] = g->p[idx];
    g->d[ghstidx] = g->d[idx];
    g->e[ghstidx] = g->e[idx];
#ifdef SEFL_GRAVITY
    g->phi[ghstidx] = g->phi[idx];
#endif
    g->v1[ghstidx] = g->v1[idx];
    g->v2[ghstidx] = g->v2[idx];
    g->v3[ghstidx] = g->v3[idx];
#ifdef MHD
    g->emf1[ghstidx] = g->emf1[idx];
    g->emf2[ghstidx] = g->emf2[idx];
    g->emf3[ghstidx] = g->emf3[idx];
#endif
  }
  // normal velocity
  /*for (int j=c->nghost; j<c->N2+c->nghost; j++)
  {
    ghstidx = j*c->full2;
    g->v1[ghstidx] = g->v1[ghstidx - c->N1];
  }*/
}

#endif //BC_P1

#ifdef BC_P2
static void periodicT(Consts* c, Grid* g)
{
  // start1 and end1 are shifted by 1 for periodic conditions
  // -- work from primary variables
  for (int j=c->N2+c->nghost; j<c->full2; j++) for (int i=c->nghost; i<c->N1+c->nghost; i++)
  {
    ghstidx = j*c->full1 + i; idx = ghstidx - c->N2*c->full1;
    g->p[ghstidx] = g->p[idx];
    g->d[ghstidx] = g->d[idx];
    g->e[ghstidx] = g->e[idx];
#ifdef SEFL_GRAVITY
    g->phi[ghstidx] = g->phi[idx];
#endif
    g->v1[ghstidx] = g->v1[idx];
    g->v2[ghstidx] = g->v2[idx];
    g->v3[ghstidx] = g->v3[idx];
#ifdef MHD
    g->emf1[ghstidx] = g->emf1[idx];
    g->emf2[ghstidx] = g->emf2[idx];
    g->emf3[ghstidx] = g->emf3[idx];
#endif
  }
  // normal velocity
  /*for (int i=c->nghost; i<c->N1+c->nghost; i++)
  {
    ghstidx = (c->full2-1)*c->full1 + i;
    g->v2[ghstidx] = g->v2[ghstidx - c->N2*c->full1];
  }*/
}

static void periodicB(Consts* c, Grid* g)
{
  // start and end variables shifted by 1 for periodic conditions
  // -- work from primary variables
  for (int j=0; j<c->nghost; j++) for (int i=c->nghost; i<c->N1+c->nghost; i++)
  {
    ghstidx = j*c->full1 + i; idx = ghstidx + c->N2*c->full1;
    g->p[ghstidx] = g->p[idx];
    g->d[ghstidx] = g->d[idx];
    g->e[ghstidx] = g->e[idx];
#ifdef SEFL_GRAVITY
    g->phi[ghstidx] = g->phi[idx];
#endif
    g->v1[ghstidx] = g->v1[idx];
    g->v2[ghstidx] = g->v2[idx];
    g->v3[ghstidx] = g->v3[idx];
#ifdef MHD
    g->emf1[ghstidx] = g->emf1[idx];
    g->emf2[ghstidx] = g->emf2[idx];
    g->emf3[ghstidx] = g->emf3[idx];
#endif
  }
  // normal velocity
  /*for (int i=c->nghost; i<c->N1+c->nghost; i++)
  {
    g->v2[i] = g->v2[i - c->N2*c->full1];
  }*/
}

#endif //BC_P2

#ifdef BC_P1
static void pL(Consts* c, double* ary)
{
  // start1 and end1 are shifted by 1 for periodic conditions
  // -- work from primary variables
  for (int j=c->nghost; j<c->N2+c->nghost; j++) for (int i=0; i<c->nghost; i++)
  {
    ghstidx = j*c->full1 + i; idx = ghstidx + c->N1;
    ary[ghstidx] = ary[idx];
  }
}

static void pR(Consts* c, double* ary)
{
  // start and end variables shifted by 1 for periodic conditions
  // -- work from primary variables
  for (int j=c->nghost; j<c->N2+c->nghost; j++) for (int i=c->N1+c->nghost; i<c->full1; i++)
  {
    ghstidx = j*c->full1 + i; idx = ghstidx - c->N1;
    ary[ghstidx] = ary[idx];
  }
}

#endif //BC_P1

#ifdef BC_P2
static void pT(Consts* c, double* ary)
{
  // start1 and end1 are shifted by 1 for periodic conditions
  // -- work from primary variables
  for (int j=c->N2+c->nghost; j<c->full2; j++) for (int i=c->nghost; i<c->N1+c->nghost; i++)
  {
    ghstidx = j*c->full1 + i; idx = ghstidx - c->N2*c->full1;
    ary[ghstidx] = ary[idx];
  }
}

static void pB(Consts* c, double* ary)
{
  // start and end variables shifted by 1 for periodic conditions
  // -- work from primary variables
  for (int j=0; j<c->nghost; j++) for (int i=c->nghost; i<c->N1+c->nghost; i++)
  {
    ghstidx = j*c->full1 + i; idx = ghstidx + c->N2*c->full1;
    ary[ghstidx] = ary[idx];
  }
}

#endif //BC_P2
