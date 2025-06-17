#include <stdio.h>

#include "../include/petsc.h"

/*
void petsc_dmda_init(DM *da);
void petsc_matrixvector_init(Mat *A, Vec *x, Vec *b, DM da,
                             double *xpt, double *ypt, double *zpt, double *XC, double *YC, double *ZC);
void petsc_ksp_init(KSP *ksp, Mat A);
void petsc_solve(double ***phi, DM da, Vec *x, Vec *b, KSP ksp,
                 double *xpt, double *ypt, double *zpt, double ***uhat, double ***vhat, double ***what);
*/

int main(int argc, char **argv) {

    // PETSc initialization
    PetscInitialize(&argc, &argv, NULL, NULL);
    
    // TODO: petsc_dmda_init -> MPIx, MPIy, MPIz?
    petsc_dmda_init(&da);
    petsc_matrixvector_init(&A, &x, &b, da, xpt, ypt, zpt, XC, YC, ZC);
    petsc_ksp_init(&ksp, A);
    petsc_solve();

    // Clean up
    VecDestroy(&x);
    VecDestroy(&b);
    DMDestroy(&da);

    // PETSc finalization
    PetscFinalize();
    return 0;
}
