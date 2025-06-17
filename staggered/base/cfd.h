#ifndef __CFD_H__
#define __CFD_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "petsc.h"
#include "config.h"

#define SQUARE(X) ((X) * (X))

double x[NX+1];
double y[NY+1];

double xc[NX+2];
double yc[NY+2];

double u_hat[NX+1][NY+2]; // x 방향 중간 속도 변화량
double v_hat[NX+2][NY+1]; // y 방향 중간 속도 변화량

// 다음 step 계산 후 값 삭제 방지
double u_temp[NX+1][NY+2];
double v_temp[NX+2][NY+1];
double p_temp[NX+2][NY+2];

double conxold[NX][NY];
double conyold[NX][NY];

// 초기 조건 설정 함수
void initialize(double u[NX+1][NY+2], double v[NX+2][NY+1], double p[NX+2][NY+2]);

void setGridPosition(double g);

void saveVelocity(double u[NX+1][NY+2], double v[NX+2][NY+1]);

// 압력 경계 조건 설정 함수
void setPressureBoundaryConditions(double p[NX+2][NY+2]);

// 중간 속도 갱신 함수
void updateIntermediateVelocity(double u[NX+1][NY+2], double v[NX+2][NY+1],
                                double p[NX+2][NY+2]);

void thomasMethod2D(int n, double a[], double b[], double c[], double d[][n],
                    double x[][n]);
void getLinearSolutionWithThomas(int dim, int numOfSolutions,
                                 double A[dim][dim],
                                 double RHS[][numOfSolutions],
                                 double x[][numOfSolutions]);

void updateIntermediateVelocity_fs(double u[NX+1][NY+2], double v[NX+2][NY+1],
                                   double p[NX+2][NY+2]);

// 속도 경계 조건 설정 함수
void setVelocityBoundaryConditions(double u[NX+1][NY+2], double v[NX+2][NY+1]);

// 경계 조건 설정 함수
void setBoundaryConditions(double u[NX+1][NY+2], double v[NX+2][NY+1],
                           double p[NX+2][NY+2]);

// 압력 필드 갱신 함수
void updatePressureField(double u[NX+1][NY+2], double v[NX+2][NY+1], double p[NX+2][NY+2]);

// 속도 필드 갱신 함수
void updateVelocityField(double u[NX+1][NY+2], double v[NX+2][NY+1], double p[NX+2][NY+2]);

// Absolute Error 구하기
double getAbsError(double u[NX+1][NY+2], double v[NX+2][NY+1]);

void printMatrix(int mtxNX, int mtxNY, double mtx[mtxNX][mtxNY]);

// 결과 출력 함수
void printResults(double u[NX+1][NY+2]);

#endif
