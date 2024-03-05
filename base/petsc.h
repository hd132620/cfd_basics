#ifndef __PETSC_H__
#define __PETSC_H__

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include <petscdm.h>
#include <petscdmda.h>
#include <petscksp.h>
#include <petscsys.h>
#include <petscvec.h>

#include "config.h"
// #include "../include/halo_exchange.h"

KSP ksp;
DM da;
Mat A;
Vec x_petsc, b;
PC pc;

extern int nx, ny, MPIx, MPIy, MPIz, g;

void petsc_dmda_init(DM *da);
void petsc_matrixvector_init(Mat *A, Vec *x, Vec *b, DM da,
                             double *xpt, double *ypt, double *zpt, double *XC, double *YC, double *ZC);
void petsc_ksp_init(KSP *ksp, Mat A);
void petsc_solve(double phi[NX][NY], DM da, Vec *x, Vec *b, KSP ksp,
                 double *xpt, double *ypt, double *zpt, double uhat[NX][NY], double vhat[NX][NY], double what[NX][NY]);

#endif