#include "constants.h"
#include "grids.h"

int idx, ghstidx;

// function pointers to boundary conditions updators
// Left, Right, Top, and Bottom
void (*bcL)() = nullptr;
void (*bcR)() = nullptr;
void (*bcT)() = nullptr;
void (*bcB)() = nullptr;

// forward declaration of private functions
#ifdef BC_P1
static void periodicL(Consts*, Grid*);
static void periodicR(Consts*, Grid*);
#endif
#ifdef BC_P1
static void periodicT(Consts*, Grid*);
static void periodicB(Consts*, Grid*);
#endif

//-Public Functions-----------------------------------------------------------//
// bcs_init(void)
//   sets appropriate pointers to private bc functions
void bcs_init(void)
{
#ifdef BC_P1
  bcL = &periodicL;
  bcR = &periodicR;
#endif
#ifdef BC_P2
  bcT = &periodicT;
  bcB = &periodicB;
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
}

//-Private Functions----------------------------------------------------------//
#ifdef BC_P1
static void periodicL(Consts* c, Grid* g)
{
  // start1 and end1 are shifted by 1 for periodic conditions
  // -- work from primary variables
  for (int j=c->nghost; j<c->N2+c->nghost; j++) for (int i=0; i<c->nghost; i++)
  {
    ghstidx = j*c->full1 + i; idx = ghostidx + c->N1;
    g->p[ghstidx] = g->p[idx];
    g->d[ghstidx] = g->d[idx];
    g->e[ghstidx] = g->e[idx];
#ifdef SEFL_GRAVITY
    g->phi[ghstidx] = g->phi[idx];
#endif
    g->v2[ghstidx] = g->v2[idx];
#ifdef MHD
    g->emf1[ghstidx] = g->emf1[idx];
    g->emf2[ghstidx] = g->emf2[idx];
    g->emf3[ghstidx] = g->emf3[idx];
#endif
  }
  // normal velocity
  for (int j=c->nghost; j<c->nghost+c->N2; j++)
  {
    ghstidx = j*c->full2;
    g->v1[ghstidx] = g->v1[ghstidx + c->N1];
  }
}

static void periodicR(Consts* c, Grid* g)
{
  // start and end variables shifted by 1 for periodic conditions
  // -- work from primary variables
  for (int j=c->nghost; j<c->N2+c->nghost; j++) for (int i=c->N1+c->nghost; i<c->full1; i++)
  {
    ghstidx = j*c->full1 + i; idx = ghostidx - c->N1;
    g->p[ghstidx] = g->p[idx];
    g->d[ghstidx] = g->d[idx];
    g->e[ghstidx] = g->e[idx];
#ifdef SEFL_GRAVITY
    g->phi[ghstidx] = g->phi[idx];
#endif
    g->v2[ghstidx] = g->v2[idx];
#ifdef MHD
    g->emf1[ghstidx] = g->emf1[idx];
    g->emf2[ghstidx] = g->emf2[idx];
    g->emf3[ghstidx] = g->emf3[idx];
#endif
  }
  // normal velocity
  for (int j=c->nghost; j<c->N2+c->nghost; j++)
  {
    ghstidx = j*c->full2;
    g->v1[ghstidx] = g->v1[ghstidx - c->N1];
  }
}

#endif //BC_P1

#ifdef BC_P2
static void periodicT(Consts* c, Grid* g)
{
  // start1 and end1 are shifted by 1 for periodic conditions
  // -- work from primary variables
  for (int j=c->N2+c->nghost; j<c->full2; j++) for (int i=c->nghost; i<c->N1+c->nghost; i++)
  {
    ghstidx = j*c->full1 + i; idx = ghostidx + c->N2*c->full1;
    g->p[ghstidx] = g->p[idx];
    g->d[ghstidx] = g->d[idx];
    g->e[ghstidx] = g->e[idx];
#ifdef SEFL_GRAVITY
    g->phi[ghstidx] = g->phi[idx];
#endif
    g->v1[ghstidx] = g->v1[idx];
#ifdef MHD
    g->emf1[ghstidx] = g->emf1[idx];
    g->emf2[ghstidx] = g->emf2[idx];
    g->emf3[ghstidx] = g->emf3[idx];
#endif
  }
  // normal velocity
  for (int i=c->nghost; i<c->N1+c->nghost; i++)
  {
    ghstidx = (c->full2-1)*c->full1 + i;
    g->v2[ghstidx] = g->v2[ghstidx - c*N2*c->full1];
  }
}

static void periodicB(Consts* c, Grid* g)
{
  // start and end variables shifted by 1 for periodic conditions
  // -- work from primary variables
  for (int j=0; j<c->nghost; j++) for (int i=c->nghost; i<c->N1+c->nghost; i++)
  {
    ghstidx = j*c->full1 + i; idx = ghostidx + c->N2*c->full1;
    g->p[ghstidx] = g->p[idx];
    g->d[ghstidx] = g->d[idx];
    g->e[ghstidx] = g->e[idx];
#ifdef SEFL_GRAVITY
    g->phi[ghstidx] = g->phi[idx];
#endif
    g->v1[ghstidx] = g->v1[idx];
#ifdef MHD
    g->emf1[ghstidx] = g->emf1[idx];
    g->emf2[ghstidx] = g->emf2[idx];
    g->emf3[ghstidx] = g->emf3[idx];
#endif
  }
  // normal velocity
  for (int i=c->nghost; i<c->N1+c->nghost; i++)
  {
    g->v2[i] = g->v2[i - c->N2*c->full1];
  }
}

#endif //BC_P2
