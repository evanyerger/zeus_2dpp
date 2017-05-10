#ifndef _SOURCE_
#define _SOURCE_

#include "constants.h"
#include "grids.h"

//double *v1na, *v2na, *v1nb, *v2nb, *q1, *q2, *enb;

void source_init(int);

void source_step(Consts* c, Grid* g);

void source_destruct(void);

#endif
