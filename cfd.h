#ifndef __CFD_H__
#define __CFD_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NX 33  // x 방향 격자점 수
#define NY 33  // y 방향 격자점 수
#define L 1.0   // 계산 영역의 길이
#define U_WALL 1.0  // 상단 벽면 속도 (lid-driven)
#define TMAX 5000   // 최대 시간 스텝 수
#define DT 0.001    // 시간 스텝 크기
#define DX (L / (NX - 1))  // x 방향 격자 간격
#define DY (L / (NY - 1))  // y 방향 격자 간격
#define Re 10.0 // Reynolds number 
#define MAX_ITER 1000  // 압력 보정 최대 반복 횟수
#define DP 0.01  // 압력 보정 강도

double u_hat[NX][NY];  // x 방향 중간 속도 변화량
double v_hat[NX][NY];  // y 방향 중간 속도 변화량

// 초기 조건 설정 함수
void initialize(double u[NX][NY], double v[NX][NY], double p[NX][NY]);

// 압력 경계 조건 설정 함수
void setPressureBoundaryConditions(double p[NX][NY]);

// 중간 속도 갱신 함수
void updateIntermediateVelocity(double u[NX][NY], double v[NX][NY], double p[NX][NY]);

// 속도 경계 조건 설정 함수
void setVelocityBoundaryConditions(double u[NX][NY], double v[NX][NY]);

// 경계 조건 설정 함수
void setBoundaryConditions(double u[NX][NY], double v[NX][NY], double p[NX][NY]);

// 압력 필드 갱신 함수
void updatePressureField(double u[NX][NY], double v[NX][NY], double p[NX][NY]);

// 속도 필드 갱신 함수
void updateVelocityField(double u[NX][NY], double v[NX][NY], double p[NX][NY]);

void printMatrix(double mtx[NX][NY]);

// 결과 출력 함수
void printResults(double u[NX][NY]);

#endif
