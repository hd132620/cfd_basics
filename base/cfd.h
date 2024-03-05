#ifndef __CFD_H__
#define __CFD_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "petsc.h"
#include "config.h"

#define SQUARE(X) ((X) * (X))

double x[NX];
double y[NY];

double xc[NX];
double yc[NY];

double u_hat[NX][NY]; // x 방향 중간 속도 변화량
double v_hat[NX][NY]; // y 방향 중간 속도 변화량

// 다음 step 계산 후 값 삭제 방지
double u_temp[NX][NY];
double v_temp[NX][NY];
double p_temp[NX][NY];

double conxold[NX][NY];
double conyold[NX][NY];

// 초기 조건 설정 함수
void initialize(double u[NX][NY], double v[NX][NY], double p[NX][NY]);

void setGridPosition(double g);

void saveVelocity(double u[NX][NY], double v[NX][NY]);

// 압력 경계 조건 설정 함수
void setPressureBoundaryConditions(double p[NX][NY]);

// 중간 속도 갱신 함수
void updateIntermediateVelocity(double u[NX][NY], double v[NX][NY],
                                double p[NX][NY]);

void thomasMethod2D(int n, double a[], double b[], double c[], double d[][n],
                    double x[][n]);
void getLinearSolutionWithThomas(int dim, int numOfSolutions,
                                 double A[dim][dim],
                                 double RHS[][numOfSolutions],
                                 double x[][numOfSolutions]);

void updateIntermediateVelocity_fs(double u[NX][NY], double v[NX][NY],
                                   double p[NX][NY]);

// 속도 경계 조건 설정 함수
void setVelocityBoundaryConditions(double u[NX][NY], double v[NX][NY]);

// 경계 조건 설정 함수
void setBoundaryConditions(double u[NX][NY], double v[NX][NY],
                           double p[NX][NY]);

// 압력 필드 갱신 함수
void updatePressureField(double u[NX][NY], double v[NX][NY], double p[NX][NY]);

// 속도 필드 갱신 함수
void updateVelocityField(double u[NX][NY], double v[NX][NY], double p[NX][NY]);

// Absolute Error 구하기
double getAbsError(double u[NX][NY], double v[NX][NY]);

void printMatrix(double mtx[NX][NY]);

// 결과 출력 함수
void printResults(double u[NX][NY]);

#endif
