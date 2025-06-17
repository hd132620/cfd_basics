#include "cfd.h"

void initialize(double u[NX+1][NY+2], double v[NX+2][NY+1], double p[NX+2][NY+2]) {
  // 초기 조건 설정
  for (int i = 0; i < NX; ++i) {
    for (int j = 0; j < NY; ++j) {
      u[i][j] = 0.0;
      v[i][j] = 0.0;
      p[i][j] = 0.0;
    }
  }
  setGridPosition(G);
  setBoundaryConditions(u, v, p); // 경계 조건 설정
}

void setGridPosition(double g) {
  for (int i = 0; i < NX + 1; ++i) {
    // x[i] = (1.0 - tanh(g * (1.0 - 2.0 * i / NX)) / tanh(g)) / 2 * L;
    x[i] = (double)i / (NX);
  }
  for (int j = 0; j < NY + 1; ++j) {
    // y[j] = (1.0 - tanh(g * (1.0 - 2.0 * j / NY)) / tanh(g)) / 2 * L;
    y[j] = (double)j / (NY);
  }
  for (int i = 0; i < NX; ++i) {
    xc[i] = (x[i] + x[i + 1]) / 2.0;
  }
  for (int j = 0; j < NY; ++j) {
    yc[j] = (y[j] + y[j + 1]) / 2.0;
  }
}

void saveVelocity(double u[NX+1][NY+2], double v[NX+2][NY+1]) {
  for (int i = 0; i < NX + 1; ++i) {
    for (int j = 0; j < NY + 2; ++j) {
      u_temp[i][j] = u[i][j];
    }
  }
  for (int i = 0; i < NX + 2; ++i) {
    for (int j = 0; j < NY + 1; ++j) {
      v_temp[i][j] = v[i][j];
    }
  }
}

// 경계 조건 설정 함수
void setPressureBoundaryConditions(double p[NX+2][NY+2]) {
  for (int i = 1; i < NX + 1; ++i) {
    // 상단 벽면: dp/dy = 0 (Neumann 경계 조건)
    p[i][NY + 1] = p[i][NY];
    // 하단 벽면: dp/dy = 0 (Neumann 경계 조건)
    p[i][0] = p[i][1];
  }
  p[0][0] = 0.0;
  p[0][NY + 1] = 0.0;

  for (int j = 1; j < NY + 1; ++j) {
    // 왼쪽 벽면: dp/dx = 0 (Neumann 경계 조건)
    p[0][j] = p[1][j];
    // 오른쪽 벽면: dp/dx = 0 (Neumann 경계 조건)
    p[NX + 1][j] = p[NX][j];
  }
  p[NX + 1][0] = 0.0;
  p[NX + 1][NY + 1] = 0.0;
}

void setVelocityBoundaryConditions(double u[NX+1][NY+2], double v[NX+2][NY+1]) {
  double gamma_top = 0.0;
  double gamma_bottom = 0.0;
  double gamma_west = 0.0;
  double gamma_east = 0.0;

  // gamma_west = (xc[1] - x[1]) / (xc[1] - xc[0]);
  // gamma_east = (x[NX - 1] - xc[NX - 2]) / (xc[NX - 1] - xc[NX - 2]);
  // gamma_top = (y[NY - 1] - yc[NY - 2]) / (yc[NY - 1] - yc[NY - 2]);
  // gamma_bottom = (yc[1] - y[1]) / (yc[1] - yc[0]);
  for (int i = 0; i < NX + 2; ++i) {
    // 상단 벽면
    v[i][NY] = 0.0;
    // 하단 벽면: u=0, v=0
    v[i][0] = 0.0;
  }
  for (int i = 0; i < NX + 1; ++i) {
    // 상단 벽면 (lid-driven): u=U_WALL, v=0
    u[i][NY + 1] = 2 * U_WALL - u[i][NY];
    // 하단
    u[i][0] = - u[i][1];
  }
  for (int j = 0; j < NY + 1; ++j) {
    // 왼쪽 벽면: v=0
    v[0][j] = - v[1][j];
    // 오른쪽
    v[NX + 1][j] = - v[NX][j];
  }
  for (int j = 0; j < NY + 2; ++j) {
    // 왼쪽 벽면: u=0
    u[0][j] = 0.0;
    // 오른쪽
    u[NX][j] = 0.0;
  }
}
// 경계 조건 설정 함수
void setBoundaryConditions(double u[NX+1][NY+2], double v[NX+2][NY+1],
                           double p[NX+2][NY+2]) {
  setPressureBoundaryConditions(p);
  setVelocityBoundaryConditions(u, v);
}

// 중간 속도 갱신 함수
void updateIntermediateVelocity(double u[NX+1][NY+2], double v[NX+2][NY+1],
                                double p[NX+2][NY+2]) {

  // 중간 속도 갱신
  for (int i = 1; i <= NX - 1; ++i) {
    for (int j = 1; j <= NY; ++j) {
      u_hat[i][j] =
          u[i][j] +
          DT * (-((SQUARE(u[i][j]) - SQUARE(u[i - 1][j])) / (1.0/(double)NX) +
                  (u[i][j] * v[i][j] - u[i][j - 1] * v[i][j - 1]) / (1.0/(double)NY)) +
                    ((u[i + 1][j] - 2.0 * u[i][j] + u[i - 1][j]) *
                         ((double)(NX*NX)) +
                     (u[i][j + 1] - 2.0 * u[i][j] + u[i][j - 1]) *
                         ((double)(NY*NY))) / Re);
      // if (u_hat[i][j] != 0.0)
      //   printf("%d %d => %lf\n",i, j, u_hat[i][j]);
      // if(isnan(u_hat[i][j]) != 0) {
      //   printf("%i %i -> %lf %lf\n", i, j, (1.0/NX),
      //               (
      //                    (yc[j + 1] - yc[j])));
      // }
    }
  }
  for (int i = 1; i <= NX; ++i) {
    for (int j = 1; j <= NY - 1; ++j) {
      v_hat[i][j] =
          v[i][j] +
          DT * (-((u[i][j] * v[i][j] - u[i - 1][j] * v[i - 1][j]) / (1.0/(double)NX) +
                  (SQUARE(v[i][j]) - SQUARE(v[i][j - 1])) / (1.0/(double)NY)) +
                1.0 / Re *
                    (
                      (v[i + 1][j] - 2.0 * v[i][j] + v[i - 1][j]) /
                         (1.0/(double)(NX*NX)) +
                     (v[i][j + 1] - 2.0 * v[i][j] + v[i][j - 1]) /
                         (1.0/(double)(NY*NY))));
    }
  }
}

void thomasMethod2D(int n, double a[], double b[], double c[], double d[][n],
                    double x[][n]) {
  // Forward Elimination
  for (int i = 1; i < n; i++) {
    double m = a[i] / b[i - 1];
    b[i] -= m * c[i - 1];
    for (int j = 0; j < n; j++) {
      c[i] -= m * c[i - 1];
    }
    for (int j = 0; j < n; j++) {
      d[i][j] -= m * d[i - 1][j];
    }
  }

  // Backward Substitution
  for (int k = 0; k < n; k++) {
    x[n - 1][k] = d[n - 1][k] / b[n - 1];
    for (int i = n - 2; i >= 0; i--) {
      double sum = d[i][k];
      for (int j = i + 1; j < n; j++) {
        sum -= c[i] * x[j][k];
      }
      x[i][k] = sum / b[i];
    }
  }
}

void getLinearSolutionWithThomas(int dim, int numOfSolutions,
                                 double A[dim][dim],
                                 double RHS[][numOfSolutions],
                                 double x[][numOfSolutions]) {
  double *a = (double *)malloc(dim * sizeof(double *));
  double *b = (double *)malloc(dim * sizeof(double *));
  double *c = (double *)malloc(dim * sizeof(double *));

  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      if (i - j == 1)
        a[i] = -x[i][j];
      else if (i == j)
        b[i] = 1 - x[i][j];
      else if (i - j == -1)
        c[i] = -x[i][j];
    }
  }
  a[0] = 0.0;
  c[dim - 1] = 0.0;
  thomasMethod2D(numOfSolutions, a, b, c, RHS, x);
  free(a);
  free(b);
  free(c);
}

void updateIntermediateVelocity_fs(double u[NX+1][NY+2], double v[NX+2][NY+1],
                                   double p[NX+2][NY+2]) {
  // convection variables
  double conx[NX][NY];
  double cony[NX][NY];

  // diffusion variables
  // double diffx[NX][NY];
  // double diffy[NX][NY];

  // A variables;
  double A1u[NX][NY];
  double A2u[NX][NY];

  double A1v[NX][NY];
  double A2v[NX][NY];

  // RHS variables with f(conv, diff)
  double rhsx[NX][NY];
  double rhsy[NX][NY];

  for (int i = 1; i < NX - 1; ++i) {
    for (int j = 1; j < NY - 1; ++j) {
      // x
      // convection term
      conx[i][j] = -DT / 2.0 *
                   (u[i][j] * (u[i + 1][j] - u[i][j]) / (xc[i + 1] - xc[i]) +
                    v[i][j] * (u[i][j + 1] - u[i][j]) / (yc[j + 1] - yc[j]));
      // A term
      A1u[i][j] =
          DT / 2.0 / Re *
          (u[i - 1][j] / (xc[i + 1] - xc[i]) / (x[i] - x[i - 1]) +
           u[i][j] / (xc[i + 1] - xc[i]) *
               ((-1.0 / (x[i + 1] - x[i])) + (-1.0 / (x[i] - x[i - 1]))) +
           u[i + 1][j] / (xc[i + 1] - xc[i]) / (x[i + 1] - x[i]));
      A2u[i][j] =
          DT / 2.0 / Re *
          (u[i - 1][j] / (yc[i + 1] - yc[i]) / (y[i] - y[i - 1]) +
           u[i][j] / (yc[i + 1] - yc[i]) *
               ((-1.0 / (y[i + 1] - y[i])) + (-1.0 / (y[i] - y[i - 1]))) +
           u[i + 1][j] / (yc[i + 1] - yc[i]) / (y[i + 1] - y[i]));
      // RHS calculation
      rhsx[i][j] =
          1.5 * conx[i][j] - 0.5 * conxold[i][j] + 2 * (A1u[i][j] + A2u[i][j]);
      conxold[i][j] = conx[i][j];
      // y
      // convection term
      cony[i][j] = -DT / 2.0 *
                   (u[i][j] * (v[i + 1][j] - v[i][j]) / (xc[i + 1] - xc[i]) +
                    v[i][j] * (v[i][j + 1] - v[i][j]) / (yc[j + 1] - yc[j]));
      // A term
      A1v[i][j] =
          DT / 2.0 / Re *
          (v[i - 1][j] / (xc[i + 1] - xc[i]) / (x[i] - x[i - 1]) +
           v[i][j] / (xc[i + 1] - xc[i]) *
               ((-1.0 / (x[i + 1] - x[i])) + (-1.0 / (x[i] - x[i - 1]))) +
           v[i + 1][j] / (xc[i + 1] - xc[i]) / (x[i + 1] - x[i]));
      A2v[i][j] =
          DT / 2.0 / Re *
          (v[i - 1][j] / (yc[i + 1] - yc[i]) / (y[i] - y[i - 1]) +
           v[i][j] / (yc[i + 1] - yc[i]) *
               ((-1.0 / (y[i + 1] - y[i])) + (-1.0 / (y[i] - y[i - 1]))) +
           v[i + 1][j] / (yc[i + 1] - yc[i]) / (y[i + 1] - y[i]));
      // RHS calculation
      rhsy[i][j] =
          1.5 * cony[i][j] - 0.5 * conyold[i][j] + 2 * (A1v[i][j] + A2v[i][j]);
      conyold[i][j] = cony[i][j];
    }
  }

  // TODO: 풀이에서 Boundary 제외
  // (1 - A1)(1 - A2)X = RHS (X = u_hat - u)
  //
  // (1 - A1)X'= RHS
  //     (NX, NX) x (NX, 1) NY개= (NX, 1) NY개
  double rhsSolutionX_[NX][NY];
  getLinearSolutionWithThomas(NX, NY, A1u, rhsx, rhsSolutionX_);

  // (1 - A2)X = X'
  //     (NX, NX) x (NX, 1) NY개 = (NX, 1) NY개
  double rhsSolutionX[NX][NY];
  getLinearSolutionWithThomas(NX, NY, A2u, rhsx, rhsSolutionX);

  // (1 - A1)(1 - A2)Y = RHS (Y = v_hat - v)
  //
  // (1 - A1)Y'= RHS
  //     (NX, NX) x (NX, 1) NY개= (NX, 1) NY개
  double rhsSolutionY_[NX][NY];
  getLinearSolutionWithThomas(NX, NY, A1v, rhsy, rhsSolutionY_);

  // (1 - A2)Y = Y'
  //     (NX, NX) x (NX, 1) NY개 = (NX, 1) NY개
  double rhsSolutionY[NX][NY];
  getLinearSolutionWithThomas(NX, NY, A2v, rhsy, rhsSolutionY);

  for (int i = 1; i < NX - 1; ++i) {
    for (int j = 1; j < NY - 1; ++j) {
      u_hat[i][j] = u[i][j] + rhsSolutionX[i][j];
      v_hat[i][j] = v[i][j] + rhsSolutionY[i][j];
    }
  }
}

int its = 0;
void updatePressureField(double u[NX+1][NY+2], double v[NX+2][NY+1], double p[NX+2][NY+2]) {
  // double dp; // 압력 보정값
  // double D, M, N;
  // double minus_M;
  // double vij_minus_vijminus1;

  // // 압력 보정을 위한 반복 수행
  // for (int it = 0; it < 10000; ++it) {
  //   // 압력 기울기 항 계산
  //   for (int i = 1; i < NX - 1; ++i) {
  //     for (int j = 1; j < NY - 1; ++j) {
  //       // vij_minus_vijminus1 = (v[i][j] - v[i][j-1]) == 0.0 ? 1.0 : (v[i][j] -
  //       // v[i][j-1]);
  //       vij_minus_vijminus1 = (v[i][j] - v[i][j - 1]) + 0.000001; // epsilon
  //       D = (u_hat[i][j] - u_hat[i - 1][j]) / (x[i] - x[i - 1]) +
  //           (v_hat[i][j] - v_hat[i][j - 1]) / (y[j] - y[j - 1]);
  //       D /= DT;
  //       minus_M =
  //           -1.0 / (xc[i + 1] - xc[i]) + (-1.0 / (xc[i] - xc[i - 1])) +
  //           (u[i][j] - u[i - 1][j]) / vij_minus_vijminus1 *
  //               (-1.0 / (yc[j + 1] - yc[j]) + (-1.0) / (yc[j] - yc[j - 1]));
  //       M = -minus_M;
  //       N = p[i + 1][j] / (xc[i + 1] - xc[i]) +
  //           p[i - 1][j] / (xc[i] - xc[i - 1]) +
  //           ((u[i][j] - u[i - 1][j]) / vij_minus_vijminus1) *
  //               (p[i][j + 1] / (yc[j + 1] - yc[j]) +
  //                p[i][j - 1] / (yc[j] - yc[j - 1]));

  //       dp = (N - (u[i][j] - u[i - 1][j]) * D) / M;

  //       // 압력 갱신
  //       // p[i][j] = dp;
  //       p[i][j] = (1 - DP) * p[i][j] + DP * dp;
  //     }
  //   }
  // }
  // its += 1;
  
  petsc_dmda_init(&da);
  petsc_matrixvector_init(&A, &x_petsc, &b, da, x, y, NULL, xc, yc, NULL);
  petsc_ksp_init(&ksp, A);
  petsc_solve(p, da, &x_petsc, &b, ksp, x, y, NULL, u_hat, v_hat, NULL);

  // Clean up
  VecDestroy(&x_petsc);
  VecDestroy(&b);
  DMDestroy(&da);
}

void updateVelocityField(double u[NX+1][NY+2], double v[NX+2][NY+1], double p[NX+2][NY+2]) {
  // 속도 갱신
  for (int i = 1; i <= NX - 1; ++i) {
    for (int j = 1; j <= NY; ++j) {
      u[i][j] =
          u_hat[i][j] - DT * (p[i][j] - p[i - 1][j]) / (1.0/(double)NX);
    }
  }
  for (int i = 1; i <= NX; ++i) {
    for (int j = 1; j <= NY - 1; ++j) {
      v[i][j] =
          v_hat[i][j] - DT * (p[i][j] - p[i][j - 1]) / (1.0/(double)NY);
    }
  }
}

double getAbsError(double u[NX+1][NY+2], double v[NX+2][NY+1]) {
  double nom = 0.0;
  double denom = 0.0;
  double error = 0.0;

  for (int i = 1; i < NX; ++i) {
    for (int j = 1; j < NY + 1; ++j) {
      nom += fabs(u[i][j] - u_temp[i][j]);
      denom += fabs(u[i][j]);
    }
  }
  for (int i = 1; i < NX + 1; ++i) {
    for (int j = 1; j < NY; ++j) {
      nom += fabs(v[i][j] - v_temp[i][j]);
      denom += fabs(v[i][j]);
    }
  }
  // printf("nom : %.4lf, denom : %.4lf\n", nom, denom);
  return nom / denom;
}

void printMatrix(int mtxNX, int mtxNY, double mtx[mtxNX][mtxNY]) {
  // 결과 출력
  int divide = 4;
  for (int j = mtxNY - 1; j >= 0; j--) {
    for (int i = 0; i < mtxNX; i++) {
      if (i % divide == 0 && j % divide == 0) {
        printf("%3.2e ", mtx[i][j]);
      }
    }
    if (j % divide == 0)
      puts("");
  }
}

void printResults(double u[NX+1][NY+2]) {
  // 결과 출력
  for (int i = 0; i < NX; i++) {
    for (int j = 0; j < NY; j++) {
      printf("%.2e ", u[i][j]);
    }
    puts("");
  }
}
