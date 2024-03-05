#include "../include/petsc.h"

// void coordinate_init(double *xpt, double *ypt, double *zpt, double *XC, double *YC, double *ZC) {
//     for (int i = 0; i < nx; ++i) {
//         xpt[i] = (1.0 - tanh(g * (1.0 - 2.0 * i / (nx - 1))) / tanh(g)) / 2 * L;
//     }
//     for (int j = 0; j < ny; ++j) {
//         ypt[j] = (1.0 - tanh(g * (1.0 - 2.0 * j / (ny - 1))) / tanh(g)) / 2 * L;
//     }
//     for (int k = 0; k < nz; ++k) {
//         zpt[k] = (1.0 - tanh(g * (1.0 - 2.0 * k / (nz - 1))) / tanh(g)) / 2 * L;
//     }
//     for (int i = 0; i < nx; ++i) {
//         XC[i] = (xpt[i] + xpt[i + 1]) / 2;
//     }
//     for (int j = 0; j < ny; ++j) {
//         YC[j] = (ypt[j] + ypt[j + 1]) / 2;
//     }
//     for (int k = 0; k < nz; ++k) {
//         ZC[k] = (zpt[k] + zpt[k + 1]) / 2;
//     }
// }

void petsc_dmda_init(DM *da)
{
    // create dmda
    PetscErrorCode ierr;

    PetscInt *lx;
    PetscInt *ly;
    PetscInt *lz;

    lx = (PetscInt *)calloc(MPIx, sizeof(PetscInt));
    ly = (PetscInt *)calloc(MPIy, sizeof(PetscInt));
    lz = (PetscInt *)calloc(MPIz, sizeof(PetscInt));

    int i;

    for (i = 0; i < MPIx; i++)
    {
        int tempsize;
        tempsize = (int)round((double)nx / MPIx);
        // if (i==(MPIx-1)) lx[i]=nx-(MPIx-1)*tempsize-1;
        if (i == (MPIx - 1))
            lx[i] = nx - (MPIx - 1) * tempsize;
        // else if (i==0) lx[i] = tempsize-1;
        else
            lx[i] = tempsize;
    }

    for (i = 0; i < MPIy; i++)
    {
        int tempsize;
        tempsize = (int)round((double)ny / MPIy);
        // if (i==(MPIy-1)) ly[i]=ny-(MPIy-1)*tempsize-1;
        if (i == (MPIy - 1))
            ly[i] = ny - (MPIy - 1) * tempsize;
        // else if (i==0) ly[i] = tempsize-1;
        else
            ly[i] = tempsize;
    }

    for (i = 0; i < MPIz; i++)
    {
        int tempsize;
        tempsize = (int)round((double)nz / MPIz);
        // if (i==(MPIz-1)) lz[i]=nz-(MPIz-1)*tempsize-1;
        if (i == (MPIz - 1))
            lz[i] = nz - (MPIz - 1) * tempsize;
        // else if (i==0) lz[i] = tempsize-1;
        else
            lz[i] = tempsize;
    }

    PetscInt dof = 1;
    PetscInt stencil_width = 1;

    ierr = DMDACreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_PERIODIC,
                        DMDA_STENCIL_STAR, nx, ny, nz, MPIx, MPIy, MPIz, dof, stencil_width, lx, ly, lz, da);
    // DMDA_STENCIL_STAR, nx-2, ny-2, nz-2, MPIx, MPIy, MPIz, dof, stencil_width, lx, ly, lz, da);

    ierr = DMSetFromOptions(*da);
    ierr = DMSetUp(*da);

    free(lx);
    free(ly);
    free(lz);
}

void petsc_matrixvector_init(Mat *A, Vec *x, Vec *b, DM da,
                             double *xpt, double *ypt, double *zpt, double *XC, double *YC, double *ZC)
{
    // create matrix, vector
    PetscErrorCode ierr;

    MatStencil matrow, matcol[7];
    PetscScalar matv[7];
    PetscInt xs, ys, zs, xm, ym, zm;

    ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);
    ierr = DMSetMatType(da, "aij");
    ierr = DMCreateMatrix(da, A);

    double ap, ae, aw, an, as, ab, af;

    // if(SHIFT_FLAG==0)
    int i, j, k;
    {
        // global index system
        for (k = 1 + zs; k < 1 + zs + zm; k++) // 1 must be added since Petsc handles only internel cells
        {
            for (j = 1 + ys; j < 1 + ys + ym; j++)
            {
                for (i = 1 + xs; i < 1 + xs + xm; i++)
                {

                    ae = (ypt[j] - ypt[j - 1]) * (zpt[k] - zpt[k - 1]) / (XC[i + 1] - XC[i]);
                    aw = (ypt[j] - ypt[j - 1]) * (zpt[k] - zpt[k - 1]) / (XC[i] - XC[i - 1]);

                    an = (xpt[i] - xpt[i - 1]) * (zpt[k] - zpt[k - 1]) / (YC[j + 1] - YC[j]);
                    as = (xpt[i] - xpt[i - 1]) * (zpt[k] - zpt[k - 1]) / (YC[j] - YC[j - 1]);

                    af = (xpt[i] - xpt[i - 1]) * (ypt[j] - ypt[j - 1]) / (ZC[k + 1] - ZC[k]);
                    ab = (xpt[i] - xpt[i - 1]) * (ypt[j] - ypt[j - 1]) / (ZC[k] - ZC[k - 1]);

                    if (i == 1)
                        aw = 0.0;
                    // if(i==nx)	ae = 0.0;

                    if (j == 1)
                        as = 0.0;
                    // if(j==ny)	an = 0.0;

                    ap = ae + aw + an + as + af + ab;

                    matrow.i = i - 1;
                    matrow.j = j - 1;
                    matrow.k = k - 1;

                    int cnt = 0;
                    matv[cnt] = ap;
                    matcol[cnt].i = i - 1;
                    matcol[cnt].j = j - 1;
                    matcol[cnt].k = k - 1; // p
                    cnt++;

                    // Neumann
                    if (i != 1)
                    {
                        matv[cnt] = -aw;
                        matcol[cnt].i = i - 1 - 1;
                        matcol[cnt].j = j - 1;
                        matcol[cnt].k = k - 1; // w
                        cnt++;
                    }

                    if (i != nx)
                    {
                        matv[cnt] = -ae;
                        matcol[cnt].i = i + 1 - 1;
                        matcol[cnt].j = j - 1;
                        matcol[cnt].k = k - 1; // e
                        cnt++;
                    }

                    // Neumann
                    if (j != 1)
                    {
                        matv[cnt] = -as;
                        matcol[cnt].i = i - 1;
                        matcol[cnt].j = j - 1 - 1;
                        matcol[cnt].k = k - 1; // s
                        cnt++;
                    }
                    // if(j!=ny-2){
                    if (j != ny)
                    {
                        matv[cnt] = -an;
                        matcol[cnt].i = i - 1;
                        matcol[cnt].j = j + 1 - 1;
                        matcol[cnt].k = k - 1; // n
                        cnt++;
                    }

                    // periodic
                    matv[cnt] = -ab;
                    matcol[cnt].i = i - 1;
                    matcol[cnt].j = j - 1;
                    matcol[cnt].k = k - 1 - 1; // b
                    cnt++;

                    matv[cnt] = -af;
                    matcol[cnt].i = i - 1;
                    matcol[cnt].j = j - 1;
                    matcol[cnt].k = k + 1 - 1; // f
                    cnt++;

                    ierr = MatSetValuesStencil(*A, 1, &matrow, cnt, matcol, matv, INSERT_VALUES);
                    // CHKERRQ(ierr);
                }
            }
        }
    }

    MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY);

    ierr = DMCreateGlobalVector(da, b);
    ierr = DMCreateGlobalVector(da, x);
}

void petsc_ksp_init(KSP *ksp, Mat A)
{
    PetscErrorCode ierr;

    ierr = KSPCreate(PETSC_COMM_WORLD, ksp);
    ierr = KSPSetOperators(*ksp, A, A);

    ierr = PCCreate(PETSC_COMM_WORLD, &pc);

    // ********** original ************
    // ierr = KSPSetType(*ksp, KSPCG);
    // ierr = KSPSetTolerances(*ksp, tol_rel, PETSC_DEFAULT, PETSC_DEFAULT, max_iter);
    // ierr = KSPSetFromOptions(*ksp);
    // ierr = KSPSetUp(*ksp);

    // ********** LU ************
    // ierr = KSPGetPC(*ksp, &pc);
    // ierr = PCSetType(pc , PCLU);
    // ierr = PCFactorSetMatSolverType(pc, MATSOLVERSUPERLU_DIST);

    // ierr = KSPSetType(*ksp, KSPCG);
    // ierr = KSPSetTolerances(*ksp, tol_rel, PETSC_DEFAULT, PETSC_DEFAULT, max_iter);
    // ierr = KSPSetFromOptions(*ksp);
    // ierr = KSPSetUp(*ksp);

    // ********** AMG ************
    ierr = KSPSetType(*ksp, KSPPREONLY);
    ierr = KSPGetPC(*ksp, &pc);
    ierr = PCSetType(pc, PCGAMG);
    ierr = PCMGSetLevels(pc, 2, NULL);

    ierr = KSPSetTolerances(*ksp, tol_rel, PETSC_DEFAULT, PETSC_DEFAULT, max_iter);
    ierr = KSPSetFromOptions(*ksp);
    ierr = KSPSetUp(*ksp);


}

void petsc_solve(double ***phi, DM da, Vec *x, Vec *b, KSP ksp,
                 double *xpt, double *ypt, double *zpt, double ***uhat, double ***vhat, double ***what)
{
    PetscErrorCode ierr;
    PetscScalar ***xx, ***bb;

    PetscInt xs, ys, zs, xm, ym, zm;
    ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);

    ierr = DMDAVecGetArray(da, *b, &bb);
    // if(myid==0) printf("check_petsc\n");

    int i, j, k;
    for (k = 1; k < zsize - 1; k++)
    {
        for (j = 1; j < ysize - 1; j++)
        {
            for (i = 1; i < xsize - 1; i++)
            {

                bb[k - 1 + zs][j - 1 + ys][i - 1 + xs] = -(uhat[i][j][k] - uhat[i - 1][j][k]) * (ypt[j + ystart] - ypt[j + ystart - 1]) * (zpt[k + zstart] - zpt[k + zstart - 1]) / dt - (vhat[i][j][k] - vhat[i][j - 1][k]) * (xpt[i + xstart] - xpt[i + xstart - 1]) * (zpt[k + zstart] - zpt[k + zstart - 1]) / dt - (what[i][j][k] - what[i][j][k - 1]) * (xpt[i + xstart] - xpt[i + xstart - 1]) * (ypt[j + ystart] - ypt[j + ystart - 1]) / dt;
            }
        }
    }
    // if(myid==0) printf("check_petsc\n");

    ierr = DMDAVecRestoreArray(da, *b, &bb);

    // if(myid==0) printf("check_petsc\n");

    VecAssemblyBegin(*b);
    VecAssemblyEnd(*b);

    // if(myid==0) printf("check_petsc\n");

    ierr = DMDAVecGetArray(da, *x, &xx);
    for (k = 1; k < zsize - 1; k++)
    {
        for (j = 1; j < ysize - 1; j++)
        {
            for (i = 1; i < xsize - 1; i++)
            {
                xx[k - 1 + zs][j - 1 + ys][i - 1 + xs] = phi[i][j][k];
            }
        }
    }

    ierr = DMDAVecRestoreArray(da, *x, &xx);
    ierr = VecAssemblyBegin(*x);
    ierr = VecAssemblyEnd(*x);

    // if(myid==0) printf("check_petsc\n");

    PetscInt its;
    PetscReal rnorm;

    ierr = KSPSetFromOptions(ksp);
    ierr = KSPSolve(ksp, *b, *x);
    ierr = KSPGetIterationNumber(ksp, &its);
    ierr = KSPGetResidualNorm(ksp, &rnorm);
    ierr = KSPGetSolution(ksp, x);

   // if (myid == 0)
     //   printf("petsc finish! %d/%e \n", its, rnorm);

    ierr = DMDAVecGetArray(da, *x, &xx); // From vector to array
    for (k = 1; k < zsize - 1; k++)
    {
        for (j = 1; j < ysize - 1; j++)
        {
            for (i = 1; i < xsize - 1; i++)
            {
                phi[i][j][k] = xx[k - 1 + zs][j - 1 + ys][i - 1 + xs];
            }
        }
    }
    ierr = DMDAVecRestoreArray(da, *x, &xx);

    double pref;
    if (myid == 0)
        pref = phi[1][1][1]; // Reference Pressure at the inlet (p_inf)
    MPI_Bcast(&pref, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (k = 0; k < zsize - 2; k++)
    {
        for (j = 0; j < ysize - 2; j++)
        {
            for (i = 0; i < xsize - 2; i++)
            {
                phi[i + 1][j + 1][k + 1] = phi[i + 1][j + 1][k + 1] - pref;
            }
        }
    }

    halo_exchange_xdir_single(phi);
    halo_exchange_ydir_single(phi);
    halo_exchange_zdir_single(phi);    

}

