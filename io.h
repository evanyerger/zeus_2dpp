#ifndef _IO_
#define _IO_
#include <string>
#include "constants.h"
#include "grids.h"
#include "timestep.h"

int check_setup(std::string);
void load_sim(std::string, Consts*, Grid*, TimeKeeper*, int);
void save_sim(std::string);
void data_dump(std::string, Consts*, Grid*, double, TimeKeeper*);
void load_variable(std::string, std::string, Consts*, double*, int);
void save_variable(std::string, std::string, Consts*, double*, int);

#endif
