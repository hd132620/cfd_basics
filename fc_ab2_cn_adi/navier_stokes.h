#ifndef NAVIER_STOKES_H
#define NAVIER_STOKES_H

#include "setting.h"

// For air properties
#define _CONST_Re 100.0  // Reynolds number
#define RHO 1.0  // kg/m^3
#define NU (1.0 / (_CONST_Re))  // m^2/s (kinematic viscosity)

extern double **u, **v, **p;

void ns_init(void);
void ns_solve(void);
void ns_finalize(void);
void ns_print_results(const char *filename_u, const char *filename_v, const char *filename_p);

#endif
