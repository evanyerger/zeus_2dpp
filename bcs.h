#ifndef _BCS_
#define _BCS_

#include "constants.h"
#include "grids.h"

void bcs_init(void);

void update_bcs(Consts*, Grid*);

void single_bc(Consts*, double*);

#endif
