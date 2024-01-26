#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cfd.h"

#define SQUARE(X) ((X) * (X))

void initialize(double u[NX][NY], double v[NX][NY], double p[NX][NY]) {
    // 초기 조건 설정
    for (int i = 1; i < NX - 1; ++i) {
        for (int j = 1; j < NY - 1; ++j) {
            u[i][j] = 0.0;
            v[i][j] = 0.0;
            p[i][j] = 0.0;
        }
    }
    setBoundaryConditions(u, v, p);  // 경계 조건 설정
}

// 경계 조건 설정 함수
void setPressureBoundaryConditions(double p[NX][NY]) {
    for (int i = 0; i < NX; ++i) {
        // 상단 벽면: dp/dy = 0 (Neumann 경계 조건)
        p[i][NY-1] = p[i][NY-2];
        // 하단 벽면: dp/dy = 0 (Neumann 경계 조건)
        p[i][0] = p[i][1];
    }

    for (int j = 0; j < NY; ++j) {
        // 왼쪽 벽면: dp/dx = 0 (Neumann 경계 조건)
        p[0][j] = p[1][j];
        // 오른쪽 벽면: dp/dx = 0 (Neumann 경계 조건)
        p[NX-1][j] = p[NX-2][j];
    }
}

void setVelocityBoundaryConditions(double u[NX][NY], double v[NX][NY]) {
    for (int i = 0; i < NX; ++i) {
        // 하단 벽면: u=0, v=0
        u[i][0] = 0.0;
        v[i][0] = 0.0;
        // 상단 벽면 (lid-driven): u=U_WALL, v=0
        u[i][NY-1] = U_WALL;
        v[i][NY-1] = 0.0;
    }

    for (int j = 0; j < NY; ++j) {
        // 왼쪽 벽면: u=0, v=0
        u[0][j] = 0.0;
        v[0][j] = 0.0;
        // 오른쪽 벽면: u=0, v=0
        u[NX-1][j] = 0.0;
        v[NX-1][j] = 0.0;
    }
}

void setBoundaryConditions(double u[NX][NY], double v[NX][NY], double p[NX][NY]) {
    setPressureBoundaryConditions(p);
    setVelocityBoundaryConditions(u, v);
}

// 중간 속도 갱신 함수
void updateIntermediateVelocity(double u[NX][NY], double v[NX][NY], double p[NX][NY]) {

    // 중간 속도 갱신
    for (int i = 1; i < NX - 1; ++i) {
        for (int j = 1; j < NY - 1; ++j) {
            // u_hat[i][j] = u[i][j] + DT * (-((p[i + 1][j] - p[i][j]) / DX) + NU * ((u[i + 1][j] - 2.0 * u[i][j] + u[i - 1][j]) / (DX * DX) + (u[i][j + 1] - 2.0 * u[i][j] + u[i][j - 1]) / (DY * DY)));
            double u_i_half_left_j = (u[i][j] + u[i-1][j]) / 2;
            double u_i_half_right_j = (u[i][j] + u[i+1][j]) / 2;
            double u_i_j_half_left = (u[i][j] + u[i][j-1]) / 2;
            double u_i_j_half_right = (u[i][j] + u[i][j+1]) / 2;

            double v_i_half_left_j = (v[i][j] + v[i-1][j]) / 2;
            double v_i_half_right_j = (v[i][j] + v[i+1][j]) / 2;
            double v_i_j_half_left = (v[i][j] + v[i][j-1]) / 2;
            double v_i_j_half_right = (v[i][j] + v[i][j+1]) / 2;

            u_hat[i][j] = u[i][j]
                        + DT * ( 
                            - 1.0 * (
                                    (SQUARE(u_i_half_right_j) - SQUARE((u_i_half_left_j))) / DX
                                + (
                                    (u_i_j_half_right * v_i_half_right_j) 
                                    - (u_i_j_half_left * v_i_half_right_j)
                                ) / DY
                            )
                            + 1.0 / Re * (
                                ((u[i+1][j] - u[i][j]) - (u[i][j] - u[i-1][j])) / DX / DX
                                + ((u[i][j+1] - u[i][j]) - (u[i][j] - u[i][j-1])) / DY / DY
                            )
                        );
            v_hat[i][j] = v[i][j]
                        + DT * ( 
                            - 1.0 * (
                                    (SQUARE(v_i_j_half_right) - SQUARE((v_i_j_half_left))) / DY
                                + (
                                    (u_i_j_half_right * v_i_half_right_j) 
                                    - (u_i_j_half_right * v_i_half_left_j)
                                ) / DX
                            )
                            + 1.0 / Re * (
                                ((v[i+1][j] - v[i][j]) - (v[i][j] - v[i-1][j])) / DX / DX
                                + ((v[i][j+1] - v[i][j]) - (v[i][j] - v[i][j-1])) / DY / DY
                            )
                        );
        }
    }
}
int its = 0;
void updatePressureField(double u[NX][NY], double v[NX][NY], double p[NX][NY]) {
    double dp; // 압력 보정값
    double D, M, N;
    double vij_minus_vijminus1;

    // 압력 보정을 위한 반복 수행
    for (int it = 0; it < MAX_ITER; ++it) {
        // 압력 기울기 항 계산
        for (int i = 1; i < NX - 1; ++i) {
            for (int j = 1; j < NY - 1; ++j) {
                vij_minus_vijminus1 = (v[i][j] - v[i][j-1]) == 0.0 ? 1.0 : (v[i][j] - v[i][j-1]);
                D = (u_hat[i][j] - u_hat[i-1][j]) / DX
                    + (v_hat[i][j] - v_hat[i][j-1]) / DY;
                D /= DT;
                M = (1.0 / DX)
                    + (1.0 / DY)
                        * (u[i][j] - u[i-1][j]) / vij_minus_vijminus1;
                N = (p[i+1][j] / DX + p[i-1][j] / DX)
                    + (p[i][j+1] / DY + p[i][j-1] / DY)
                        * ((u[i][j] - u[i-1][j]) / vij_minus_vijminus1);
                    
                dp = N - (u[i][j] - u[i-1][j]) * D / M;

                // 압력 갱신
                p[i][j] = (1 - DP) * p[i][j] + DP * dp;
            }
        }
    }
    its += 1;
}

void updateVelocityField(double u[NX][NY], double v[NX][NY], double p[NX][NY]) {
    // 속도 갱신
    for (int i = 1; i < NX - 1; ++i) {
        for (int j = 1; j < NY - 1; ++j) {
            u[i][j] = u_hat[i][j] - DT * (p[i + 1][j] - p[i - 1][j]) / (2.0 * DX);
            v[i][j] = v_hat[i][j] - DT * (p[i][j + 1] - p[i][j - 1]) / (2.0 * DY);
        }
    }
}

void printMatrix(double mtx[NX][NY]) {
    // 결과 출력
    for (int j = NY - 1; j >= 0; j--) {
        for (int i = 0; i < NX; i++) {
            printf("%.1lf ", mtx[i][j]);
        }
        puts("");
    }
}

void printResults(double u[NX][NY]) {
    // 결과 출력
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            printf("%.1lf ", u[i][j]);
        }
        puts("");
    }
}
