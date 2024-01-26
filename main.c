#include <stdio.h>
#include "cfd.h"

int main(void) {
    double u[NX][NY];  // x 방향 속도
    double v[NX][NY];  // y 방향 속도
    double p[NX][NY];  // 압력

    puts("실행");

    initialize(u, v, p);  // 초기 조건 설정
    puts("초기화 완료");
    puts("Step 시작");
    for (int t = 0; t < MAX_ITER; ++t) {
        printf("Step %d : %lfs\n", t + 1, DT * t);
        updateIntermediateVelocity(u, v, p);  // 중간 속도 갱신
        updatePressureField(u, v, p);    // 압력 갱신
        updateVelocityField(u, v, p);    // 속도 갱신
        setBoundaryConditions(u, v, p);  // 경계 조건 설정
        printMatrix(u);
    }
    return 0;
}
