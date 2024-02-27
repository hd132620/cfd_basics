#include "cfd.h"
#include <stdio.h>

int main(void) {
  double u[NX][NY]; // x 방향 속도
  double v[NX][NY]; // y 방향 속도
  double p[NX][NY]; // 압력
  const double tol = 1.0;
  double absError = 0;
  int t = 0;

  puts("실행");

  initialize(u, v, p); // 초기 조건 설정
  puts("초기화 완료");
  puts("Step 시작");
  do {
    puts("------------------------");
    printf("Step %d : %lfs\n", t + 1, DT * t);
    // printMatrix(u);
    saveVelocity(u, v);
    updateIntermediateVelocity_fs(u, v, p); // 중간 속도 갱신
    updatePressureField(u, v, p);           // 압력 갱신
    updateVelocityField(u, v, p);           // 속도 갱신
    setBoundaryConditions(u, v, p);         // 경계 조건 설정
    printMatrix(v);
    absError = getAbsError(u, v);
    printf("Tol : %.2lf%% < %.2lf%%\n ", absError * 100, tol * 100);
    t += 1;
  } while (absError <= tol && t <= TMAX);
  return 0;
}
