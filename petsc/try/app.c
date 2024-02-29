#include <stdlib.h>
#include <petscdmda.h>
#include <petscksp.h>
#include <petscsys.h>
#include <petscvec.h>
#include "../../include/config.h"

#define DT 0.001
#define NX 32
#define NY 32
#define MAX_ITER 100
#define DP 0.9
#define Re 1000
#define U_WALL 1.0

PetscReal x[NX], y[NY], xc[NX], yc[NY];
double conxold[NY][NX], conyold[NY][NX];
PC pc;

void printMatrix(double mtx[NY][NX]);

// Function declarations
void updateIntermediateVelocity_fs(DM da, Vec u_vec, Vec v_vec,
                                             Vec p_vec, Vec u_hat_vec,
                                             Vec v_hat_vec);
void updatePressureField(DM da, Vec u_vec, Vec v_vec, Vec p_vec);
void updateVelocityField(DM da, Vec u_vec, Vec v_vec, Vec p_vec, Vec u_hat_vec,
                         Vec v_hat_vec);

void setPressureBoundaryConditions(DM da, Vec p_vec);
void setVelocityBoundaryConditions(DM da, Vec u_vec, Vec v_vec);
void setBoundaryConditions(DM da, Vec u_vec, Vec v_vec, Vec p_vec);

PetscErrorCode petsc_matrixvector_init(Mat *A, Vec *x, Vec *b, DM da,
                              double A1_fs[][NX], double A2_fs[][NX]);
PetscErrorCode petsc_ksp_init(KSP *ksp, Mat A);
PetscErrorCode petsc_solve(double solution[][NX], DM da, Vec *x, Vec *b, KSP ksp, double rhs[][NX]);

void petsc_dmda_init(DM *da);

int main(int argc, char **argv) {
  PetscInt Nx = NX, Ny = NY;
  PetscInt g = 2;
  PetscReal Lx = 1.0, Ly = 1.0;

  for (int i = 0; i < Nx; ++i) {
    x[i] = (1.0 - tanh(g * (1.0 - 2.0 * i / (Nx - 1))) / tanh(g)) / 2 * Lx;
  }
  for (int j = 0; j < Ny; ++j) {
    y[j] = (1.0 - tanh(g * (1.0 - 2.0 * j / (Ny - 1))) / tanh(g)) / 2 * Ly;
  }
  for (int i = 0; i < Nx; ++i) {
    xc[i] = (x[i] + x[i + 1]) / 2;
  }
  for (int j = 0; j < Ny; ++j) {
    yc[j] = (y[j] + y[j + 1]) / 2;
  }

  // PETSc initialization
  PetscInitialize(&argc, &argv, NULL, NULL);

  // Create DMDA
  DM da;
  petsc_dmda_init(&da);

  // Create Petsc vectors for u, v, and p
  Vec u_vec, v_vec, p_vec, u_hat_vec, v_hat_vec;
  DMCreateGlobalVector(da, &u_vec);
  DMCreateGlobalVector(da, &v_vec);
  DMCreateGlobalVector(da, &p_vec);
  DMCreateGlobalVector(da, &u_hat_vec);
  DMCreateGlobalVector(da, &v_hat_vec);

  double temp[NY][NX];

  // loop
  for (int iter = 0; iter < MAX_ITER; iter++) {
    setBoundaryConditions(da, u_vec, v_vec, p_vec);
    // DMDAVecGetArray(da, u_vec, temp);
    // printMatrix(temp);
    // updateIntermediateVelocity_fs(da, u_vec, v_vec, p_vec, u_hat_vec,
    //                               v_hat_vec);
    // updatePressureField(da, u_vec, v_vec, p_vec);
    // updateVelocityField(da, u_vec, v_vec, p_vec, u_hat_vec, v_hat_vec);
  }

  // Clean up
  VecDestroy(&u_vec);
  VecDestroy(&v_vec);
  VecDestroy(&p_vec);
  DMDestroy(&da);

  // PETSc finalization
  PetscFinalize();
  return 0;
}

void printMatrix(double mtx[NY][NX]) {
  // 결과 출력
  int divide = 1;
  for (int j = NY - 1; j >= 0; j--) {
    for (int i = 0; i < NX; i++) {
      if (i % divide == 0 && j % divide == 0) {
        printf("%.1lf ", mtx[i][j]);
      }
    }
    if (j % divide == 0)
      puts("");
  }
}

void updateIntermediateVelocity_fs(DM da, Vec u_vec, Vec v_vec,
                                             Vec p_vec, Vec u_hat_vec,
                                             Vec v_hat_vec) {
  
  KSP ksp;
  PetscReal **u, **v, **p, **u_hat, **v_hat;
  DMDALocalInfo info;
  PetscInt i, j, mx, my;

  DMDAVecGetArray(da, u_vec, &u);
  DMDAVecGetArray(da, v_vec, &v);
  DMDAVecGetArray(da, p_vec, &p);
  DMDAVecGetArray(da, u_hat_vec, &u_hat);
  DMDAVecGetArray(da, v_hat_vec, &v_hat);

  DMDAGetInfo(da, NULL, &mx, &my, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
              NULL, NULL, NULL);

  DMDAGetLocalInfo(da, &info);

  // convection variables
  double conx[NY][NX];
  double cony[NY][NX];

  // A variables;
  double A1u[NY][NX];
  double A2u[NY][NX];

  double A1v[NY][NX];
  double A2v[NY][NX];

  // RHS variables with f(conv, diff)
  double rhsx[NY][NX];
  double rhsy[NY][NX];

  for (j = info.ys; j < info.ys + info.ym; j++) {
    for (i = info.xs; i < info.xs + info.xm; i++) {
      // x
      conx[j][i] = -DT / 2.0 *
                   (u[j][i] * (u[j][i + 1] - u[j][i]) / (xc[i + 1] - xc[i]) +
                    v[j][i] * (u[j + 1][i] - u[j][i]) / (yc[j + 1] - yc[j]));
      A1u[j][i] =
          DT / 2.0 / Re *
          (u[j][i - 1] / (xc[i + 1] - xc[i]) / (x[i] - x[i - 1]) +
           u[j][i] / (xc[i + 1] - xc[i]) *
               ((-1.0 / (x[i + 1] - x[i])) + (-1.0 / (x[i] - x[i - 1]))) +
           u[j][i + 1] / (xc[i + 1] - xc[i]) / (x[i + 1] - x[i]));
      A2u[j][i] =
          DT / 2.0 / Re *
          (u[j - 1][i] / (yc[i + 1] - yc[i]) / (y[j] - y[j - 1]) +
           u[j][i] / (yc[i + 1] - yc[i]) *
               ((-1.0 / (y[j + 1] - y[j])) + (-1.0 / (y[j] - y[j - 1]))) +
           u[j + 1][i] / (yc[i + 1] - yc[i]) / (y[j + 1] - y[j]));
      rhsx[j][i] =
          1.5 * conx[j][i] - 0.5 * conxold[j][i] + 2 * (A1u[j][i] + A2u[j][i]);
      conxold[j][i] = conx[j][i];
      // y
      cony[j][i] = -DT / 2.0 *
                   (u[j][i] * (v[j + 1][i] - v[j][i]) / (xc[i + 1] - xc[i]) +
                    v[j][i] * (v[j][i + 1] - v[j][i]) / (yc[j + 1] - yc[j]));
      A1v[j][i] =
          DT / 2.0 / Re *
          (v[j - 1][i] / (xc[i + 1] - xc[i]) / (x[i] - x[i - 1]) +
           v[j][i] / (xc[i + 1] - xc[i]) *
               ((-1.0 / (x[i + 1] - x[i])) + (-1.0 / (x[i] - x[i - 1]))) +
           v[j + 1][i] / (xc[i + 1] - xc[i]) / (x[i + 1] - x[i]));
      A2v[j][i] =
          DT / 2.0 / Re *
          (v[j - 1][i] / (yc[i + 1] - yc[i]) / (y[j] - y[j - 1]) +
           v[j][i] / (yc[i + 1] - yc[i]) *
               ((-1.0 / (y[j + 1] - y[j])) + (-1.0 / (y[j] - y[j - 1]))) +
           v[j + 1][i] / (yc[i + 1] - yc[i]) / (y[j + 1] - y[j]));
      rhsy[j][i] =
          1.5 * cony[j][i] - 0.5 * conyold[j][i] + 2 * (A1v[j][i] + A2v[j][i]);
      conyold[j][i] = cony[j][i];
    }
  }

  Mat A;
  Vec x, b;
  double rhsSolutionX[NY][NX], rhsSolutionY[NY][NX];

  petsc_matrixvector_init(&A, &x, &b, da, A1u, A2u);
  petsc_ksp_init(&ksp, A);
  petsc_solve(rhsSolutionX, da, &x, &b, ksp, rhsx);

  petsc_matrixvector_init(&A, &x, &b, da, A1v, A2v);
  petsc_ksp_init(&ksp, A);
  petsc_solve(rhsSolutionY, da, &x, &b, ksp, rhsy);

  // double _rhsx[NY][NX];
  // double _rhsy[NY][NX];

  // petsc_matrixvector_init(&A, &x, &b, A2u);
  // petsc_ksp_init(&ksp, A);
  // petsc_solve();

  // petsc_matrixvector_init(&A, &x, &b, A1v);
  // petsc_ksp_init(&ksp, A);
  // petsc_solve();

  // petsc_matrixvector_init(&A, &x, &b, A2v);
  // petsc_ksp_init(&ksp, A);
  // petsc_solve();

  // (1 - A1)(1 - A2)X = RHS (X = u_hat - u)

  // // (1 - A2)X = X'

  // // (1 - A1)(1 - A2)Y = RHS (Y = v_hat - v)

  // // (1 - A2)Y = Y'

  for (int i = 1; i < NX - 1; ++i) {
    for (int j = 1; j < NY - 1; ++j) {
      u_hat[j][i] = u[j][i] + rhsSolutionX[j][i];
      v_hat[j][i] = v[j][i] + rhsSolutionY[j][i];
    }
  }

  // Restore arrays
  DMDAVecRestoreArray(da, u_vec, &u);
  DMDAVecRestoreArray(da, v_vec, &v);
  DMDAVecRestoreArray(da, p_vec, &p);
  DMDAVecRestoreArray(da, u_hat_vec, &u_hat);
  DMDAVecRestoreArray(da, v_hat_vec, &v_hat);

}

void updatePressureField(DM da, Vec u_vec, Vec v_vec, Vec p_vec) {
  PetscInt i, j, Ma, Na, mstart, nstart, mend, nend;

  DMDAGetInfo(da, NULL, &Ma, &Na, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
              NULL, NULL, NULL);

  DMDAGetCorners(da, &mstart, &nstart, NULL, &mend, &nend, NULL);

  // Local arrays for velocity and pressure
  PetscScalar **uLocal, **vLocal, **pLocal;
  DMDAVecGetArray(da, u_vec, &uLocal);
  DMDAVecGetArray(da, v_vec, &vLocal);
  DMDAVecGetArray(da, p_vec, &pLocal);

  double dp; // Pressure correction value
  double D, M, N;
  double minus_M;
  double vij_minus_vijminus1;

  // Perform iterations for pressure correction
  for (int it = 0; it < MAX_ITER; ++it) {
    // Calculate pressure gradient term
    for (j = nstart; j < nend; ++j) {
      for (i = mstart; i < mend; ++i) {
        // vij_minus_vijminus1 = (vLocal[j][i] - vLocal[j][i-1]) == 0.0 ? 1.0 :
        // (vLocal[j][i] - vLocal[j][i-1]);
        vij_minus_vijminus1 =
            (vLocal[j][i] - vLocal[j][i - 1]) + 0.000001; // epsilon
        D = (uLocal[j][i] - uLocal[j][i - 1]) / (x[i] - x[i - 1]) +
            (vLocal[j][i] - vLocal[j - 1][i]) / (y[j] - y[j - 1]);
        D /= DT;
        minus_M =
            -1.0 / (xc[i + 1] - xc[i]) + (-1.0 / (xc[i] - xc[i - 1])) +
            (uLocal[j][i] - uLocal[j][i - 1]) / vij_minus_vijminus1 *
                (-1.0 / (yc[j + 1] - yc[j]) + (-1.0) / (yc[j] - yc[j - 1]));
        M = -minus_M;
        N = pLocal[j + 1][i] / (xc[i + 1] - xc[i]) +
            pLocal[j - 1][i] / (xc[i] - xc[i - 1]) +
            ((vLocal[j][i] - vLocal[j - 1][i]) / vij_minus_vijminus1) *
                (pLocal[j][i + 1] / (yc[j + 1] - yc[j]) +
                 pLocal[j][i - 1] / (yc[j] - yc[j - 1]));

        dp = (N - (vLocal[j][i] - vLocal[j - 1][i]) * D) / M;

        // Update pressure
        pLocal[j][i] = (1 - DP) * pLocal[j][i] + DP * dp;
      }
    }
  }

  // Restore arrays
  // DMDAVecRestoreArray(da, u_vec, &uLocal);
  // DMDAVecRestoreArray(da, v_vec, &vLocal);
  DMDAVecRestoreArray(da, p_vec, &pLocal);
}

void updateVelocityField(DM da, Vec u_vec, Vec v_vec, Vec p_vec, Vec u_hat_vec,
                         Vec v_hat_vec) {

  DMDALocalInfo info;
  PetscInt mx, my;
  PetscReal **u, **v, **p, **u_hat, **v_hat;

  DMDAGetLocalInfo(da, &info);
  DMDAGetInfo(da, NULL, &mx, &my, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
              NULL, NULL, NULL);

  DMDAVecGetArray(da, u_vec, &u);
  DMDAVecGetArray(da, v_vec, &v);
  DMDAVecGetArray(da, p_vec, &p);
  DMDAVecGetArray(da, u_hat_vec, &u_hat);
  DMDAVecGetArray(da, v_hat_vec, &v_hat);

  // Loop through local elements
  for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {
    for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
      u[j][i] =
          u_hat[j][i] - DT * (p[j][i] - p[j][i - 1]) / (xc[i] - xc[i - 1]);
      v[j][i] =
          v_hat[j][i] - DT * (p[j][i] - p[j - 1][i]) / (yc[j] - yc[j - 1]);
    }
  }

  DMDAVecRestoreArray(da, u_vec, &u);
  DMDAVecRestoreArray(da, v_vec, &v);
}

// 경계 조건 설정 함수
void setPressureBoundaryConditions(DM da, Vec p_vec) {
  DMDALocalInfo info;

  DMDAGetLocalInfo(da, &info);

  PetscScalar **pLocal;
  DMDAVecGetArray(da, p_vec, &pLocal);

  // 상단 벽면: dp/dy = 0 (Neumann 경계 조건)
  for (int i = info.xs; i < info.xs + info.xm; i++) {
    pLocal[info.ys + info.ym - 1][i] = pLocal[info.ys + info.ym - 2][i];
  }

  // 하단 벽면: dp/dy = 0 (Neumann 경계 조건)
  for (int i = info.xs; i < info.xs + info.xm; i++) {
    pLocal[info.ys][i] = pLocal[info.ys + 1][i];
  }

  // 왼쪽 벽면: dp/dx = 0 (Neumann 경계 조건)
  for (int j = info.ys; j < info.ys + info.ym; j++) {
    pLocal[j][info.xs] = pLocal[j][info.xs + 1];
  }

  // 오른쪽 벽면: dp/dx = 0 (Neumann 경계 조건)
  for (int j = info.ys; j < info.ys + info.ym; j++) {
    pLocal[j][info.xs + info.xm - 1] = pLocal[j][info.xs + info.xm - 2];
  }

  // Restore the array
  DMDAVecRestoreArray(da, p_vec, &pLocal);
}

void setVelocityBoundaryConditions(DM da, Vec u_vec, Vec v_vec) {

  DMDALocalInfo info;
  DMDAGetLocalInfo(da, &info);

  PetscScalar **uArray, **vArray;
  DMDAVecGetArray(da, u_vec, &uArray);
  DMDAVecGetArray(da, v_vec, &vArray);

  PetscReal gamma_top, gamma_bottom, gamma_west, gamma_east;

  // Set velocity boundary conditions
  for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
    // Top boundary (lid-driven): u=U_WALL, v=0
    gamma_top = (vArray[info.ys + info.ym - 2][i] - yc[NY - 2]) /
                (yc[NY - 2] - yc[NY - 3]);
    uArray[info.ys + info.ym - 1][i] =
        (U_WALL - (1 - gamma_top) * uArray[info.ys + info.ym - 2][i]) /
        gamma_top;
    vArray[info.ys + info.ym - 2][i] = 0.0;
    // Bottom boundary: u=0, v=0
    gamma_bottom = (vArray[info.ys][i] - yc[0]) / (yc[1] - yc[0]);
    uArray[info.ys] =
        -gamma_bottom / (1 - gamma_bottom) * uArray[info.ys + 1][i];
    vArray[info.ys][i] = 0.0;
  }

  for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {
    // Left boundary: u=0, v=0
    gamma_west = (uArray[j][info.xs] - xc[0]) / (xc[1] - xc[0]);
    uArray[j][info.xs] = 0.0;
    vArray[j][info.xs] =
        -gamma_west / (1 - gamma_west) * vArray[j][info.xs + 1];
    // Right boundary: u=0, v=0
    gamma_east = (uArray[j][info.xs + info.xm - 2] - xc[NX - 2]) /
                 (xc[NX - 2] - xc[NX - 3]);
    uArray[j][info.xs + info.xm - 2] = 0.0;
    vArray[j][info.xs + info.xm - 1] =
        -(1 - gamma_east) / gamma_east * vArray[j][info.xs + info.xm - 2];
  }

  // Restore arrays
  DMDAVecRestoreArray(da, u_vec, &uArray);
  DMDAVecRestoreArray(da, v_vec, &vArray);
}

void setBoundaryConditions(DM da, Vec u_vec, Vec v_vec, Vec p_vec) {
  setPressureBoundaryConditions(da, p_vec);
  setVelocityBoundaryConditions(da, u_vec, v_vec);
}


PetscErrorCode petsc_matrixvector_init(Mat *A, Vec *x, Vec *b, DM da, double A1_fs[][NX], double A2_fs[][NX])
{
    // create matrix, vector
    PetscErrorCode ierr;

    MatStencil matrow, matcol[7];
    PetscScalar matv[7];
    PetscInt xs, ys, zs, xm, ym, zm;

    ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);
    ierr = DMSetMatType(da, "aij");
    ierr = DMCreateMatrix(da, A);

    // if(SHIFT_FLAG==0)
    int i, j;
    {
        // global index system
            for (j = ys; j < ys + ym; j++)
            {
                for (i = xs; i < xs + xm; i++)
                {
                    matrow.i = i;
                    matrow.j = j;

                    int cnt = 0;

                    if (i == j) {
                      matv[cnt] = 1 - A1_fs[j][i] - A2_fs[j][i];
                      matcol[cnt].i = i;
                      matcol[cnt].j = j;
                      cnt++;
                    } else if (abs(i - j) == 1) {
                      matv[cnt] = -A1_fs[j][i] - A2_fs[j][i];
                      matcol[cnt].i = i;
                      matcol[cnt].j = j;
                      cnt++;
                    }

                    ierr = MatSetValuesStencil(*A, 1, &matrow, cnt, matcol, matv, INSERT_VALUES);
                    // CHKERRQ(ierr);
                }
            }
        
    }

    MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY);

    ierr = DMCreateGlobalVector(da, b);
    ierr = DMCreateGlobalVector(da, x);
    CHKERRQ(ierr);
    return ierr;
}

PetscErrorCode petsc_ksp_init(KSP *ksp, Mat A)
{
    PetscErrorCode ierr;
    PetscReal tol_rel = 0.0001;

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

    ierr = KSPSetTolerances(*ksp, tol_rel, PETSC_DEFAULT, PETSC_DEFAULT, MAX_ITER);
    ierr = KSPSetFromOptions(*ksp);
    ierr = KSPSetUp(*ksp);
    CHKERRQ(ierr);
    return ierr;
}

PetscErrorCode petsc_solve(double solution[][NX], DM da, Vec *x, Vec *b, KSP ksp, double rhs[][NX])
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

            bb[j - 1 + ys][i - 1 + xs] = rhs[j][i];
        }
    }
    // if(myid==0) printf("check_petsc\n");

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
                xx[j - 1 + ys][i - 1 + xs] = solution[i][j];
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
        for (j = 1; j < NY - 1; j++)
        {
            for (i = 1; i < NX - 1; i++)
            {
                solution[i][j] = xx[j - 1 + ys][i - 1 + xs];
            }
        }
    ierr = DMDAVecRestoreArray(da, *x, &xx);
    CHKERRQ(ierr); 
    return ierr;
}


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