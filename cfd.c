#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "cfd.h"

#define SQUARE(X) ((X) * (X))

void initialize(double u[NX][NY], double v[NX][NY], double p[NX][NY]) {
    // 초기 조건 설정
    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            u[i][j] = 0.0;
            v[i][j] = 0.0;
            p[i][j] = 0.0;
        }
    }
    setGridPosition(G);
    setBoundaryConditions(u, v, p);  // 경계 조건 설정
}

void setGridPosition(double g) {
    for (int i = 0; i < NX; ++i) {
        x[i] = (1.0 - tanh(g * (1.0 - 2.0 * i / (NX - 1))) / tanh(g)) / 2 * L;
    }
    for (int j = 0; j < NY; ++j) {
        y[j] = (1.0 - tanh(g * (1.0 - 2.0 * j / (NY - 1))) / tanh(g)) / 2 * L;
    }
}

void saveVelocity(double u[NX][NY], double v[NX][NY]) {
    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            u_temp[i][j] = u[i][j];
            v_temp[i][j] = v[i][j];
        }
    }
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
    double gamma_top = 0.0;
    double gamma_bottom = 0.0;
    double gamma_west = 0.0;
    double gamma_east = 0.0;

    for (int i = 0; i < NX; ++i) {
        // 상단 벽면 (lid-driven): u=U_WALL, v=0
        gamma_top = (v[i][NY-2] - y[NY-2]) / (y[NY-1] - y[NY-2]);
        u[i][NY-1] = (U_WALL - (1 - gamma_top) * u[i][NY-2]) / gamma_top;
        v[i][NY-2] = 0.0;
        // 하단 벽면: u=0, v=0
        gamma_bottom = (v[i][0] - y[0]) / (y[1] - y[0]);
        u[i][0] = -gamma_bottom / (1 - gamma_bottom) * u[i][1];
        v[i][0] = 0.0;
    }


    for (int j = 0; j < NY; ++j) {
        // 왼쪽 벽면: u=0, v=0
        gamma_west = (u[0][j] - x[0]) / (x[1] - x[0]);
        u[0][j] = 0.0;
        v[0][j] = -gamma_west / (1 - gamma_west) * v[1][j];
        // 오른쪽 벽면: u=0, v=0
        gamma_east = (u[NX-2][j] - x[NX-2]) / (x[NX-1] - x[NX-2]);
        u[NX-2][j] = 0.0;
        v[NX-1][j] = -(1 - gamma_east) / gamma_east * v[NX-2][j];
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
            u_hat[i][j] = u[i][j]
                        + DT * ( 
                            - (
                                u[i][j] * (u[i+1][j] - u[i][j]) / (x[i+1] - x[i])
                                + v[i][j] * (u[i][j+1] - u[i][j]) / (y[j+1] - y[j])
                            )
                            + 1.0 / Re * (
                                ((u[i+1][j] - u[i][j]) / (x[i+1] - x[i]) - (u[i][j] - u[i-1][j]) / (x[i] - x[i-1])) / (x[i] - x[i-1])
                                + ((u[i][j+1] - u[i][j]) / (y[j+1] - y[j]) - (u[i][j] - u[i][j-1]) / (y[j] - y[j-1])) / (y[j] - y[j-1])
                            )
                        );
            v_hat[i][j] = v[i][j]
                        + DT * ( 
                            - (
                                u[i][j] * (v[i+1][j] - v[i][j]) / (x[i+1] - x[i])
                                + v[i][j] * (v[i][j+1] - v[i][j]) / (y[j+1] - y[j])
                            )
                            + 1.0 / Re * (
                                ((v[i+1][j] - v[i][j]) / (x[i+1] - x[i]) - (v[i][j] - v[i-1][j]) / (x[i] - x[i-1])) / (x[i] - x[i-1])
                                + ((v[i][j+1] - v[i][j]) / (y[j+1] - y[j]) - (v[i][j] - v[i][j-1]) / (y[j] - y[j-1])) / (y[j] - y[j-1])
                            )
                        );
        }
    }
}

void updateIntermediateVelocity_fs(double u[NX][NY], double v[NX][NY], double p[NX][NY]) {

}

int its = 0;
void updatePressureField(double u[NX][NY], double v[NX][NY], double p[NX][NY]) {
    double dp; // 압력 보정값
    double D, M, N;
    double minus_M;
    double vij_minus_vijminus1;

    // 압력 보정을 위한 반복 수행
    for (int it = 0; it < MAX_ITER; ++it) {
        // 압력 기울기 항 계산
        for (int i = 1; i < NX - 1; ++i) {
            for (int j = 1; j < NY - 1; ++j) {
                vij_minus_vijminus1 = (v[i][j] - v[i][j-1]) == 0.0 ? 1.0 : (v[i][j] - v[i][j-1]);
                D = (u_hat[i][j] - u_hat[i-1][j]) / (x[i] - x[i-1])
                    + (v_hat[i][j] - v_hat[i][j-1]) / (y[j] - y[j-1]);
                D /= DT;
                minus_M = -1.0 / (x[i+1] - x[i])
                          + (-1.0 / (x[i] - x[i-1]))
                          + (u[i][j] - u[i-1][j]) / vij_minus_vijminus1
                            * (-1.0 / (y[j+1] - y[j]) + (-1.0) / (y[j] - y[j-1]));
                M = -minus_M;
                N = p[i+1][j] / (x[i+1] - x[i])
                    + p[i-1][j] / (x[i] - x[i-1])
                    + ((u[i][j] - u[i-1][j]) / vij_minus_vijminus1)
                        * (p[i][j+1] / (y[j+1] - y[j]) + p[i][j-1] / (y[j] - y[j-1]));
                
                dp = (N - (u[i][j] - u[i-1][j]) * D) / M;

                // 압력 갱신
                // p[i][j] = dp;
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
            u[i][j] = u_hat[i][j] - DT * (p[i][j] - p[i - 1][j]) / (x[i] - x[i-1]);
            v[i][j] = v_hat[i][j] - DT * (p[i][j] - p[i][j - 1]) / (y[j] - y[j-1]);
        }
    }
}

double getAbsError(double u[NX][NY], double v[NX][NY]) {
    double nom = 0.0;
    double denom = 0.0;

    for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
            nom += fabs(u[i][j] - u_temp[i][j]);
            nom += fabs(v[i][j] - v_temp[i][j]);
            denom += fabs(u[i][j]);
            denom += fabs(v[i][j]);
        }
    }
    // printf("nom : %.4lf, denom : %.4lf\n", nom, denom);
    return nom / denom;
}

void printMatrix(double mtx[NX][NY]) {
    // 결과 출력
    int divide = 1;
    for (int j = NY - 1; j >= 0; j--) {
        for (int i = 0; i < NX; i++) {
            if (i % divide == 0 && j % divide == 0) {
                printf("%.1lf ", mtx[i][j]);
            }    
        }
        if (j % divide == 0) puts("");
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
