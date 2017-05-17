#ifndef _SOURCE_
#define _SOURCE_
#include "constants.h"
#include "grids.h"
;
void source_init(int);

void source_step(Consts*, Grid*, double);

void source_destruct(void);

#endif
