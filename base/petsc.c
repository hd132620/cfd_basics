#include "petsc.h"

int nx = NX, ny = NY;
int MPIx = 1, MPIy = 1, MPIz = 1;
int g = 2.0;

void petsc_dmda_init(DM *da)
{
    // create dmda
    PetscErrorCode ierr;

    PetscInt *lx;
    PetscInt *ly;

    lx = (PetscInt *)calloc(MPIx, sizeof(PetscInt));
    ly = (PetscInt *)calloc(MPIy, sizeof(PetscInt));

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

    PetscInt dof = 1;
    PetscInt stencil_width = 1;

    ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                        DMDA_STENCIL_STAR, nx, ny, MPIx, MPIy, dof, stencil_width, lx, ly, da);
    // DMDA_STENCIL_STAR, nx-2, ny-2, nz-2, MPIx, MPIy, MPIz, dof, stencil_width, lx, ly, lz, da);

    ierr = DMSetFromOptions(*da);
    ierr = DMSetUp(*da);
    free(lx);
    free(ly);
}

void petsc_matrixvector_init(Mat *A, Vec *x, Vec *b, DM da,
                             double *xpt, double *ypt, double *zpt, double *XC, double *YC, double *ZC)
{
    PetscErrorCode ierr;

    // create matrix, vector
    MatStencil matrow, matcol[7];
    PetscScalar matv[7];
    PetscInt xs, ys, zs, xm, ym, zm;

    ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);
    ierr = DMSetMatType(da, "aij");
    ierr = DMCreateMatrix(da, A);

    double ap, ae, aw, an, as;

    // if(SHIFT_FLAG==0)
    int i, j;
        // global index system
            for (j = 1 + ys; j < 1 + ys + ym; j++)
            {
                for (i = 1 + xs; i < 1 + xs + xm; i++)
                {
                    ae = (ypt[j] - ypt[j - 1]) / (XC[i + 1] - XC[i]);
                    aw = (ypt[j] - ypt[j - 1]) / (XC[i] - XC[i - 1]);

                    an = (xpt[i] - xpt[i - 1]) / (YC[j + 1] - YC[j]);
                    as = (xpt[i] - xpt[i - 1]) / (YC[j] - YC[j - 1]);

                    if (i == 1)
                        aw = 0.0;
                    if(i==nx)	ae = 0.0;

                    if (j == 1)
                        as = 0.0;
                    if(j==ny)	an = 0.0;

                    ap = ae + aw + an + as;

                    matrow.i = i - 1;
                    matrow.j = j - 1;

                    int cnt = 0;
                    matv[cnt] = ap;
                    matcol[cnt].i = i - 1;
                    matcol[cnt].j = j - 1;  // p
                    cnt++;

                    // Neumann
                    if (i != 1)
                    {
                        matv[cnt] = -aw;
                        matcol[cnt].i = i - 1 - 1;
                        matcol[cnt].j = j - 1;  // w
                        cnt++;
                    }

                    if (i != nx)
                    {
                        matv[cnt] = -ae;
                        matcol[cnt].i = i + 1 - 1;
                        matcol[cnt].j = j - 1; // e
                        cnt++;
                    }

                    // Neumann
                    if (j != 1)
                    {
                        matv[cnt] = -as;
                        matcol[cnt].i = i - 1;
                        matcol[cnt].j = j - 1 - 1;  // s
                        cnt++;
                    }
                    // if(j!=ny-2){
                    if (j != ny)
                    {
                        matv[cnt] = -an;
                        matcol[cnt].i = i - 1;
                        matcol[cnt].j = j + 1 - 1;  // n
                        cnt++;
                    }

                    ierr = MatSetValuesStencil(*A, 1, &matrow, cnt, matcol, matv, INSERT_VALUES);
                    // CHKERRQ(ierr);
                    //  printf("(%d, %d) ", i, j);
                }
                // printf("\n");
            }

    MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY);
    // PetscPrintf(PETSC_COMM_WORLD, "행렬 값:\n");
    // MatView(*A, PETSC_VIEWER_STDOUT_WORLD);
    // PetscInt m, n;
    // MatGetSize(*A, &m, &n);
    // PetscPrintf(PETSC_COMM_WORLD, "행렬의 행 수: %d, 열 수: %d\n", m, n);

    ierr = DMCreateGlobalVector(da, b);
    ierr = DMCreateGlobalVector(da, x);
}

void petsc_ksp_init(KSP *ksp, Mat A)
{
    PetscErrorCode ierr;

    PetscReal tol_rel = 0.0001;
    PetscInt max_iter = 10000;

    ierr = KSPCreate(PETSC_COMM_WORLD, ksp);
    ierr = KSPSetOperators(*ksp, A, A);

    ierr = PCCreate(PETSC_COMM_WORLD, &pc);

    // ********** original ************
    ierr = KSPSetType(*ksp, KSPCG);
    ierr = KSPSetTolerances(*ksp, tol_rel, PETSC_DEFAULT, PETSC_DEFAULT, max_iter);
    ierr = KSPSetFromOptions(*ksp);
    ierr = KSPSetUp(*ksp);

    // ********** LU ************
    // ierr = KSPGetPC(*ksp, &pc);
    // ierr = PCSetType(pc , PCLU);
    // ierr = PCFactorSetMatSolverType(pc, MATSOLVERSUPERLU_DIST);

    // ierr = KSPSetType(*ksp, KSPCG);
    // ierr = KSPSetTolerances(*ksp, tol_rel, PETSC_DEFAULT, PETSC_DEFAULT, max_iter);
    // ierr = KSPSetFromOptions(*ksp);
    // ierr = KSPSetUp(*ksp);

    // ********** AMG ************
    // ierr = KSPSetType(*ksp, KSPPREONLY);
    // ierr = KSPGetPC(*ksp, &pc);
    // ierr = PCSetType(pc, PCGAMG);
    // ierr = PCMGSetLevels(pc, 2, NULL);

    // ierr = KSPSetTolerances(*ksp, tol_rel, PETSC_DEFAULT, PETSC_DEFAULT, max_iter);
    // ierr = KSPSetFromOptions(*ksp);
    // ierr = KSPSetUp(*ksp);

}

void petsc_solve(double phi[NX][NY], DM da, Vec *x, Vec *b, KSP ksp,
                 double *xpt, double *ypt, double *zpt, double uhat[NX][NY], double vhat[NX][NY], double what[NX][NY])
{
    PetscErrorCode ierr;
    PetscScalar **xx, **bb;

    PetscInt xs, ys, zs, xm, ym, zm;
    ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);

    ierr = DMDAVecGetArray(da, *b, &bb);
    // if(myid==0) printf("check_petsc\n");

    int i, j;
    for (j = 1; j < NY - 1; j++)
    {
        for (i = 1; i < NX - 1; i++)
        {
            bb[j - 1 + ys][i - 1 + xs] = -(uhat[i][j] - uhat[i - 1][j]) * (ypt[j + ys] - ypt[j + ys - 1]) / DT - (vhat[i][j] - vhat[i][j - 1]) * (xpt[i + xs] - xpt[i + xs - 1]) / DT;
            // printf("%.4e ", bb[j - 1 + ys][i - 1 + xs]);
        }
        // puts("");
    }
    // if(myid==0) printf("check_petsc\n");
    // printMatrix(bb);

    ierr = DMDAVecRestoreArray(da, *b, &bb);

    // if(myid==0) printf("check_petsc\n");

    VecAssemblyBegin(*b);
    VecAssemblyEnd(*b);

    // if(myid==0) printf("check_petsc\n");

    ierr = DMDAVecGetArray(da, *x, &xx);
    for (j = 1; j < NY - 1; j++)
    {
        for (i = 1; i < NX - 1; i++)
        {
            xx[j - 1 + ys][i - 1 + xs] = phi[i][j];
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

    puts("== Poisson PETSc KSP stage ==");
    printf("Iter: %d, Rnorm : %.2f\n", its, rnorm);
    puts("=============================");
   // if (myid == 0)
     //   printf("petsc finish! %d/%e \n", its, rnorm);

    ierr = DMDAVecGetArray(da, *x, &xx); // From vector to array
        for (j = 1; j < NY - 1; j++)
        {
            for (i = 1; i < NX - 1; i++)
            {
                phi[i][j] = xx[j - 1 + ys][i - 1 + xs];
            }
        }
    ierr = DMDAVecRestoreArray(da, *x, &xx);
}

