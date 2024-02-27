#include <petscdmda.h>
#include <petscksp.h>

#define DT 0.001
#define NX 32
#define NY 32
#define MAX_ITER 100
#define DP 0.9
#define Re 1000
#define U_WALL 1.0


PetscReal x[NX], y[NY], xc[NX], yc[NY];
double    conxold[NX][NY], conyold[NX][NY];


// Function declarations
PetscErrorCode  updateIntermediateVelocity_fs(DM da, Vec u_vec, Vec v_vec, Vec p_vec, Vec u_hat_vec, Vec v_hat_vec);
PetscErrorCode  updatePressureField(DM da, Vec u_vec, Vec v_vec, Vec p_vec);
void            updateVelocityField(DM da, Vec u_vec, Vec v_vec, Vec p_vec, Vec u_hat_vec, Vec v_hat_vec);

void            thomasMethod2D(int n, double a[], double b[], double c[], double d[][n], double x[][n]);
void            getLinearSolutionWithThomas(int dim, int numOfSolutions, double A[dim][dim], double RHS[][numOfSolutions], double x[][numOfSolutions]);
void            setPressureBoundaryConditions(DM da, Vec p_vec);
void            setVelocityBoundaryConditions(DM da, Vec u_vec, Vec v_vec);
void            setBoundaryConditions(DM da, Vec u_vec, Vec v_vec, Vec p_vec);


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
    xc[i] = (x[i] + x[i+1]) / 2;
  }
  for (int j = 0; j < Ny; ++j) {
    yc[j] = (y[j] + y[j+1]) / 2;
  }

  // PETSc initialization
  PetscInitialize(&argc, &argv, NULL, NULL);

  // Create DMDA
  DM da;
  DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
            DMDA_STENCIL_STAR, Nx, Ny, PETSC_DECIDE, PETSC_DECIDE,
            1, 1, NULL, NULL, &da);

  // Set up DMDA
  DMSetFromOptions(da);
  DMSetUp(da);

  // Create Petsc vectors for u, v, and p
  Vec u_vec, v_vec, p_vec, u_hat_vec, v_hat_vec;
  DMCreateGlobalVector(da, &u_vec);
  DMCreateGlobalVector(da, &v_vec);
  DMCreateGlobalVector(da, &p_vec);
  DMCreateGlobalVector(da, &u_hat_vec);
  DMCreateGlobalVector(da, &v_hat_vec);

  // loop
  for (int iter = 0; iter < MAX_ITER; iter++) {
    setBoundaryConditions(da, u_vec, v_vec, p_vec);
    updateIntermediateVelocity_fs(da, u_vec, v_vec, p_vec, u_hat_vec, v_hat_vec);
    updatePressureField(da, u_vec, v_vec, p_vec);
    updateVelocityField(da, u_vec, v_vec, p_vec, u_hat_vec, v_hat_vec);
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

PetscErrorCode updateIntermediateVelocity_fs(DM da, Vec u_vec, Vec v_vec, Vec p_vec, Vec u_hat_vec, Vec v_hat_vec) {
  PetscErrorCode ierr;
  
  PetscReal    **u, **v, **p, **u_hat, **v_hat;
  DMDALocalInfo  info;
  PetscInt     i, j, mx, my;
  PetscReal    **conx, **cony, **A1u, **A2u, **conyold;
  PetscReal    **A1v, **A2v, **rhsx, **rhsy;

  
  DMDAVecGetArray(da, u_vec, &u);
  DMDAVecGetArray(da, v_vec, &v);
  DMDAVecGetArray(da, p_vec, &p);
  DMDAVecGetArray(da, u_hat_vec, &u_hat);
  DMDAVecGetArray(da, v_hat_vec, &v_hat);

  DMDAGetInfo(da, NULL, &mx, &my, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

  DMDAGetLocalInfo(da, &info);

  conx = (double **)malloc(NX * sizeof(double *));
  cony = (double **)malloc(NX * sizeof(double *));
  A1u = (double **)malloc(NX * sizeof(double *));
  A2u = (double **)malloc(NX * sizeof(double *));
  conyold = (double **)malloc(NX * sizeof(double *));
  A1v = (double **)malloc(NX * sizeof(double *));
  A2v = (double **)malloc(NX * sizeof(double *));
  rhsx = (double **)malloc(NX * sizeof(double *));
  rhsy = (double **)malloc(NX * sizeof(double *));

  for (int i = 0; i < NX; i++) {
    conx[i] = (double *)malloc(NY * sizeof(double));
    cony[i] = (double *)malloc(NY * sizeof(double));
    A1u[i] = (double *)malloc(NY * sizeof(double));
    A2u[i] = (double *)malloc(NY * sizeof(double));
    conyold[i] = (double *)malloc(NY * sizeof(double));
    A1v[i] = (double *)malloc(NY * sizeof(double));
    A2v[i] = (double *)malloc(NY * sizeof(double));
    rhsx[i] = (double *)malloc(NY * sizeof(double));
    rhsy[i] = (double *)malloc(NY * sizeof(double));
  }

  for (j = info.ys; j < info.ys + info.ym; j++) {
    for (i = info.xs; i < info.xs + info.xm; i++) {
      // x
      conx[j][i] = -DT / 2.0 * (
        u[j][i] * (u[j][i + 1] - u[j][i]) / (xc[i + 1] - xc[i]) +
        v[j][i] * (u[j + 1][i] - u[j][i]) / (yc[j + 1] - yc[j])
      );
      A1u[j][i] = DT / 2.0 / Re * (
        u[j][i - 1] / (xc[i + 1] - xc[i]) / (x[i] - x[i - 1]) +
        u[j][i] / (xc[i + 1] - xc[i]) *
          ((-1.0 / (x[i + 1] - x[i])) + (-1.0 / (x[i] - x[i - 1]))) +
        u[j][i + 1] / (xc[i + 1] - xc[i]) / (x[i + 1] - x[i])
      );
      A2u[j][i] = DT / 2.0 / Re * (
        u[j - 1][i] / (yc[i + 1] - yc[i]) / (y[j] - y[j - 1]) +
        u[j][i] / (yc[i + 1] - yc[i]) *
          ((-1.0 / (y[j + 1] - y[j])) + (-1.0 / (y[j] - y[j - 1]))) +
        u[j + 1][i] / (yc[i + 1] - yc[i]) / (y[j + 1] - y[j])
      );
      rhsx[j][i] = 1.5 * conx[j][i] - 0.5 * conxold[j][i] + 2 * (A1u[j][i] + A2u[j][i]);
      conxold[j][i] = conx[j][i];
      // y
      cony[j][i] = -DT / 2.0 * (
        u[j][i] * (v[j + 1][i] - v[j][i]) / (xc[i + 1] - xc[i]) +
        v[j][i] * (v[j][i + 1] - v[j][i]) / (yc[j + 1] - yc[j])
      );
      A1v[j][i] = DT / 2.0 / Re * (
        v[j - 1][i] / (xc[i + 1] - xc[i]) / (x[i] - x[i - 1]) +
        v[j][i] / (xc[i + 1] - xc[i]) *
          ((-1.0 / (x[i + 1] - x[i])) + (-1.0 / (x[i] - x[i - 1]))) +
        v[j + 1][i] / (xc[i + 1] - xc[i]) / (x[i + 1] - x[i])
      );
      A2v[j][i] = DT / 2.0 / Re * (
        v[j - 1][i] / (yc[i + 1] - yc[i]) / (y[j] - y[j - 1]) +
        v[j][i] / (yc[i + 1] - yc[i]) *
          ((-1.0 / (y[j + 1] - y[j])) + (-1.0 / (y[j] - y[j - 1]))) +
        v[j + 1][i] / (yc[i + 1] - yc[i]) / (y[j + 1] - y[j])
      );
      rhsy[j][i] = 1.5 * cony[j][i] - 0.5 * conyold[j][i] + 2 * (A1v[j][i] + A2v[j][i]);
      conyold[j][i] = cony[j][i];
    }
  }

  // (1 - A1)(1 - A2)X = RHS (X = u_hat - u)
  double rhsSolutionX_[NY][NX];
  getLinearSolutionWithThomas(NY, NX, A1u, rhsx, rhsSolutionX_);

  // (1 - A2)X = X'
  double rhsSolutionX[NX][NY];
  getLinearSolutionWithThomas(NX, NY, A2u, rhsSolutionX_, rhsSolutionX);

  // (1 - A1)(1 - A2)Y = RHS (Y = v_hat - v)
  double rhsSolutionY_[NX][NY];
  getLinearSolutionWithThomas(NX, NY, A1v, rhsy, rhsSolutionY_);

  // (1 - A2)Y = Y'
  double rhsSolutionY[NX][NY];
  getLinearSolutionWithThomas(NX, NY, A2v, rhsSolutionY_, rhsSolutionY);

  for (int i = 1; i < NX - 1; ++i) {
    for (int j = 1; j < NY - 1; ++j) {
      u_hat[i][j] = u[i][j] + rhsSolutionX[i][j];
      v_hat[i][j] = v[i][j] + rhsSolutionY[i][j];
    }
  }

  // Restore arrays
  DMDAVecRestoreArray(da, u_vec, &u);
  DMDAVecRestoreArray(da, v_vec, &v);
  DMDAVecRestoreArray(da, p_vec, &p);
  DMDAVecRestoreArray(da, u_hat_vec, &u_hat);
  DMDAVecRestoreArray(da, v_hat_vec, &v_hat);

  // Free the allocated memory
  for (int i = 0; i < NX; i++) {
    free(conx[i]);
    free(cony[i]);
    free(A1u[i]);
    free(A2u[i]);
    free(conyold[i]);
    free(A1v[i]);
    free(A2v[i]);
    free(rhsx[i]);
    free(rhsy[i]);
    free(rhsSolutionX_[i]);
    free(rhsSolutionX[i]);
    free(rhsSolutionY_[i]);
    free(rhsSolutionY[i]);
  }

  free(conx);
  free(cony);
  free(A1u);
  free(A2u);
  free(conyold);
  free(A1v);
  free(A2v);
  free(rhsx);
  free(rhsy);
  free(rhsSolutionX_);
  free(rhsSolutionX);
  free(rhsSolutionY_);
  free(rhsSolutionY);

  return ierr;
}

PetscErrorCode updatePressureField(DM da, Vec u_vec, Vec v_vec, Vec p_vec) {
  PetscErrorCode ierr;
  PetscInt       i, j, Ma, Na, mstart, nstart, mend, nend;

  DMDAGetInfo(da, NULL, &Ma, &Na, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

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
        // vij_minus_vijminus1 = (vLocal[j][i] - vLocal[j][i-1]) == 0.0 ? 1.0 : (vLocal[j][i] - vLocal[j][i-1]);
        vij_minus_vijminus1 = (vLocal[j][i] - vLocal[j][i-1]) + 0.000001; // epsilon
        D = (uLocal[j][i] - uLocal[j][i-1]) / (x[i] - x[i-1])
          + (vLocal[j][i] - vLocal[j-1][i]) / (y[j] - y[j-1]);
        D /= DT;
        minus_M = -1.0 / (xc[i+1] - xc[i])
              + (-1.0 / (xc[i] - xc[i-1]))
              + (uLocal[j][i] - uLocal[j][i-1]) / vij_minus_vijminus1
              * (-1.0 / (yc[j+1] - yc[j]) + (-1.0) / (yc[j] - yc[j-1]));
        M = -minus_M;
        N = pLocal[j+1][i] / (xc[i+1] - xc[i])
          + pLocal[j-1][i] / (xc[i] - xc[i-1])
          + ((vLocal[j][i] - vLocal[j-1][i]) / vij_minus_vijminus1)
            * (pLocal[j][i+1] / (yc[j+1] - yc[j]) + pLocal[j][i-1] / (yc[j] - yc[j-1]));

        dp = (N - (vLocal[j][i] - vLocal[j-1][i]) * D) / M;

        // Update pressure
        pLocal[j][i] = (1 - DP) * pLocal[j][i] + DP * dp;
      }
    }
  }

  // Restore arrays
  // DMDAVecRestoreArray(da, u_vec, &uLocal);
  // DMDAVecRestoreArray(da, v_vec, &vLocal);
  DMDAVecRestoreArray(da, p_vec, &pLocal);

  return ierr;
}


void updateVelocityField(DM da, Vec u_vec, Vec v_vec, Vec p_vec, Vec u_hat_vec, Vec v_hat_vec) {
  
  DMDALocalInfo  info;
  PetscInt mx, my;
  PetscReal **u, **v, **p, **u_hat, **v_hat;

  DMDAGetLocalInfo(da, &info);
  DMDAGetInfo(da, NULL, &mx, &my, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

  DMDAVecGetArray(da, u_vec, &u);
  DMDAVecGetArray(da, v_vec, &v);
  DMDAVecGetArray(da, p_vec, &p);
  DMDAVecGetArray(da, u_hat_vec, &u_hat);
  DMDAVecGetArray(da, v_hat_vec, &v_hat);
  
  // Loop through local elements
  for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {
    for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
      u[j][i] = u_hat[j][i] - DT * (p[j][i] - p[j][i - 1]) / (xc[i] - xc[i-1]);
      v[j][i] = v_hat[j][i] - DT * (p[j][i] - p[j - 1][i]) / (yc[j] - yc[j-1]);
    }
  }

  DMDAVecRestoreArray(da, u_vec, &u);
  DMDAVecRestoreArray(da, v_vec, &v);
}

void thomasMethod2D(int n, double a[], double b[], double c[], double d[][n], double x[][n]) {
  // Forward Elimination
  for (int i = 1; i < n; i++) {
    double m = a[i] / b[i-1];
    b[i] -= m * c[i-1];
    for (int j = 0; j < n; j++) {
      c[i] -= m * c[i-1];
    }
    for (int j = 0; j < n; j++) {
      d[i][j] -= m * d[i-1][j];
    }
  }

  // Backward Substitution
  for (int k = 0; k < n; k++) {
    x[n-1][k] = d[n-1][k] / b[n-1];
    for (int i = n-2; i >= 0; i--) {
      double sum = d[i][k];
      for (int j = i+1; j < n; j++) {
        sum -= c[i] * x[j][k];
      }
      x[i][k] = sum / b[i];
    }
  }
}

void getLinearSolutionWithThomas(int dim, int numOfSolutions, double A[dim][dim], double RHS[][numOfSolutions], double x[][numOfSolutions]) {
  double * a = (double *)malloc(dim * sizeof(double *));
  double * b = (double *)malloc(dim * sizeof(double *));
  double * c = (double *)malloc(dim * sizeof(double *));

  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      if (i - j == 1) a[i] = -x[i][j];
      else if (i == j) b[i] = 1 - x[i][j];
      else if (i - j == -1) c[i] = -x[i][j];
    }
  }
  a[0] = 0.0; c[dim-1] = 0.0;
  thomasMethod2D(numOfSolutions, a, b, c, RHS, x);
  free(a);
  free(b);
  free(c);
}

// 경계 조건 설정 함수
void setPressureBoundaryConditions(DM da, Vec p_vec) {
  PetscErrorCode ierr;
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
  ierr = DMDAVecRestoreArray(da, p_vec, &pLocal);
  CHKERRABORT(PETSC_COMM_WORLD, ierr);
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
    gamma_top = (vArray[i][info.ys + info.ym - 2] - yc[NY-2]) / (yc[NY-2] - yc[NY-3]);
    uArray[i][info.ys + info.ym - 1] = (U_WALL - (1 - gamma_top) * uArray[i][info.ys + info.ym - 2]) / gamma_top;
    vArray[i][info.ys + info.ym - 2] = 0.0;
    // Bottom boundary: u=0, v=0
    gamma_bottom = (vArray[i][info.ys] - yc[0]) / (yc[1] - yc[0]);
    uArray[i][info.ys] = -gamma_bottom / (1 - gamma_bottom) * uArray[i][info.ys + 1];
    vArray[i][info.ys] = 0.0;
  }

  for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {
    // Left boundary: u=0, v=0
    gamma_west = (uArray[info.xs][j] - xc[0]) / (xc[1] - xc[0]);
    uArray[info.xs][j] = 0.0;
    vArray[info.xs][j] = -gamma_west / (1 - gamma_west) * vArray[info.xs + 1][j];
    // Right boundary: u=0, v=0
    gamma_east = (uArray[info.xs + info.xm - 2][j] - xc[NX-2]) / (xc[NX-2] - xc[NX-3]);
    uArray[info.xs + info.xm - 2][j] = 0.0;
    vArray[info.xs + info.xm - 1][j] = -(1 - gamma_east) / gamma_east * vArray[info.xs + info.xm - 2][j];
  }

  // Restore arrays
  DMDAVecRestoreArray(da, u_vec, &uArray);
  DMDAVecRestoreArray(da, v_vec, &vArray);
}


void setBoundaryConditions(DM da, Vec u_vec, Vec v_vec, Vec p_vec) {
  setPressureBoundaryConditions(da, p_vec);
  setVelocityBoundaryConditions(da, u_vec, v_vec);
}
