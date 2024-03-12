#pragma once
#ifndef __CONFIG_H__
#define __CONFIG_H__

#define NX 192      // x 방향 격자점 수
#define NY 192      // y 방향 격자점 수
#define L 1.0      // 계산 영역의 길이
#define U_WALL 1.0 // 상단 벽면 속도 (lid-driven)
#define TMAX 1000  // 최대 시간 스텝 수
#define DT 0.0001   // 시간 스텝 크기
#define Re 1000.0      // Reynolds number
// #define MAX_ITER 45000 // 압력 보정 최대 반복 횟수
#define DP 0.01       // 압력 보정 강도
#define G 1.0

#endif
