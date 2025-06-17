#ifndef POISSON_H
#define POISSON_H

#include <stdlib.h>
#include <math.h>
#include "setting.h"
#include "navier_stokes.h"

void poisson_pressure_solver(
    double **p,
    double **rhs,
    int     nx,
    int     ny,
    double  dx,
    double  dy,
    int     iter,
    double  tol
);

void cn_diffusion_adi(
    double **u_old, double **v_old,
    double **N_u,   double **N_u_old,
    double **N_v,   double **N_v_old,
    int nx, int ny, double dx, double dy,
    double dt_, double nu_,
    double **u_star, double **v_star
);


#endif // POISSON_H
