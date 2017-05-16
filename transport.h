#ifndef _TRANSPORT_
#define _TRANSPORT_
#include "constants.h"
#include "grids.h"

void transport_init(int);

void transport_step(Consts*, Grid*, double);

void transport_destruct(void);

#endif
