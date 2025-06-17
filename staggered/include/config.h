#ifndef __CONFIG_H__
#define __CONFIG_H__

#define NX 33      // x 방향 격자점 수
#define NY 33      // y 방향 격자점 수
#define L 1.0      // 계산 영역의 길이
#define U_WALL 1.0 // 상단 벽면 속도 (lid-driven)
#define TMAX 5000  // 최대 시간 스텝 수
#define DT 0.001   // 시간 스텝 크기
// #define DX (L / (NX - 1))  // x 방향 격자 간격
// #define DY (L / (NY - 1))  // y 방향 격자 간격
#define Re 1000.0      // Reynolds number
#define MAX_ITER 45000 // 압력 보정 최대 반복 횟수
#define DP 0.001       // 압력 보정 강도
#define G 2.0

int nx = NX, ny = NY;
int MPIx = 4, MPIy = 4, MPIz = 4;
int g = 2.0;

#endif