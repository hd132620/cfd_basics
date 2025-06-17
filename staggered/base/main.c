#include <stdio.h>
#include <petscsys.h>
#include <mpi.h>

#include "config.h"
#include "cfd.h"

int main(int argc, char *argv[]) {
  double u[NX+1][NY+2]; // x 방향 속도
  double v[NX+2][NY+1]; // y 방향 속도
  double p[NX+2][NY+2]; // 압력
  const double tol = 1.0;
  double absError = 0;
  int t = 0;

  // MPI_Init(&argc, &argv);

  PetscInitialize(&argc, &argv, NULL, NULL);

  puts("실행");

  initialize(u, v, p); // 초기 조건 설정
  puts("초기화 완료");
  puts("Step 시작");
  do {
    puts("------------------------");
    printf("Step %d : %lfs\n", t + 1, DT * t);
    saveVelocity(u, v);
    // setVelocityBoundaryConditions(u_hat, v_hat);
    updateIntermediateVelocity(u, v, p);    // 중간 속도 갱신
    updatePressureField(u, v, p);           // 압력 갱신
    // setPressureBoundaryConditions(p);
    updateVelocityField(u, v, p);           // 속도 갱신
    setBoundaryConditions(u, v, p);         // 경계 조건 설정
    absError = getAbsError(u, v);
    printf("Tol : %.6lf%% < %.2lf%%\n ", absError * 100, tol * 100);
    t += 1;
    // return 0;
  } while (absError <= tol && t <= TMAX); 
    printMatrix(NX+1, NY+2, u);

  // PETSc finalization
  PetscFinalize();
  // MPI_Finalize();
  return 0;
}
