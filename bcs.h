#include "constants.h"
#include "grids.h"

#ifndef _BCS_
#define _BCS_

void bcs_init(void);

void update_bcs(Consts*, Grid*);

#endif
