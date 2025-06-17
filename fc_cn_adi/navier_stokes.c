#include "navier_stokes.h"
#include "poisson.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// 전역 변수: u, v, p (cell-centered)
double **u = NULL, **v = NULL, **p = NULL;

void ns_init(void) {
    u = (double **)malloc(NX * sizeof(double *));
    v = (double **)malloc(NX * sizeof(double *));
    p = (double **)malloc(NX * sizeof(double *));
    for (int i = 0; i < NX; i++) {
        u[i] = (double *)calloc(NY, sizeof(double));
        v[i] = (double *)calloc(NY, sizeof(double));
        p[i] = (double *)calloc(NY, sizeof(double));
    }
    // 내부 전체 u=0, v=0, p=0
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            u[i][j] = 0.0;
            v[i][j] = 0.0;
            p[i][j] = 1.0;
        }
    }
    // ghost 셀 경계조건 (no-slip / moving lid)
    for (int j = 0; j < NY; j++) {
        u[0][j]    = -u[1][j];
        u[NX-1][j] = -u[NX-2][j];
        v[0][j]    = -v[1][j];
        v[NX-1][j] = -v[NX-2][j];
    }
    for (int i = 0; i < NX; i++) {
        u[i][0]    = -u[i][1];
        v[i][0]    = -v[i][1];
    }
    for (int i = 0; i < NX; i++) {
        u[i][NY-1] =  2.0*1.0 - u[i][NY-2];
        v[i][NY-1] = -v[i][NY-2];
    }
}

void ns_finalize(void) {
    for (int i = 0; i < NX; i++) {
        free(u[i]); free(v[i]); free(p[i]);
    }
    free(u); free(v); free(p);
}


void ns_print_results(const char *filename_u, const char *filename_v, const char *filename_p) {
    FILE *fu = fopen(filename_u, "w");
    FILE *fv = fopen(filename_v, "w");
    FILE *fp = fopen(filename_p, "w");
    double dx = 1.0/(NX-1), dy = 1.0/(NY-1);
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            double x = i*dx, y = j*dy;
            fprintf(fu, "%d %d %e %e %e\n", i, j, x, y, u[i][j]);
            fprintf(fv, "%d %d %e %e %e\n", i, j, x, y, v[i][j]);
            fprintf(fp, "%d %d %e %e %e\n", i, j, x, y, p[i][j]);
        }
        fprintf(fu, "\n"); fprintf(fv, "\n"); fprintf(fp, "\n");
    }
    fclose(fu); fclose(fv); fclose(fp);
}

void ns_solve(void) {
    double dx = 1.0/(NX-1), dy = 1.0/(NY-1);

    double **u_star   = (double **)malloc(NX * sizeof(double *));
    double **v_star   = (double **)malloc(NX * sizeof(double *));
    double **p_prime  = (double **)malloc(NX * sizeof(double *));
    double **rhs      = (double **)malloc(NX * sizeof(double *));
    double **N_u      = (double **)malloc(NX * sizeof(double *));
    double **N_u_old  = (double **)malloc(NX * sizeof(double *));
    double **N_v      = (double **)malloc(NX * sizeof(double *));
    double **N_v_old  = (double **)malloc(NX * sizeof(double *));
    double **u_old    = (double **)malloc(NX * sizeof(double *));
    for (int i = 0; i < NX; i++) {
        u_star[i]   = (double *)calloc(NY, sizeof(double));
        v_star[i]   = (double *)calloc(NY, sizeof(double));
        p_prime[i]  = (double *)calloc(NY, sizeof(double));
        rhs[i]      = (double *)calloc(NY, sizeof(double));
        N_u[i]      = (double *)calloc(NY, sizeof(double));
        N_u_old[i]  = (double *)calloc(NY, sizeof(double));
        N_v[i]      = (double *)calloc(NY, sizeof(double));
        N_v_old[i]  = (double *)calloc(NY, sizeof(double));
        u_old[i]    = (double *)calloc(NY, sizeof(double));
    }
    
    // u_fc v_fc
    double **u_fc = (double **)malloc((NX - 1) * sizeof(double *));
    for (int i = 0; i < NX - 1; i++) {
        u_fc[i] = (double *)calloc(NY, sizeof(double));
    }
    double **v_fc = (double **)malloc(NX * sizeof(double *));
    for (int i = 0; i < NX; i++) {
        v_fc[i] = (double *)calloc(NY - 1, sizeof(double));
    }

    for (int n = 1; n <= NT; n++) {
        //-----------------------------------------------------
        // 0) 이전 시점 속도 보관
        //-----------------------------------------------------
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                u_old[i][j] = u[i][j];
            }
        }

        //-----------------------------------------------------
        // 1) 대류항 계산: 2차 중심 차분 + AB2
        //-----------------------------------------------------
        for (int i = 1; i < NX-1; i++) {
            for (int j = 1; j < NY-1; j++) {
                // u-momentum convection: u * du/dx + v * du/dy
                double u_val = u[i][j];
                double v_val = v[i][j];
                double dudx = (u[i+1][j] - u[i-1][j]) / (2*dx);
                double dudy = (u[i][j+1] - u[i][j-1]) / (2*dy);
                N_u[i][j] = u_val * dudx + v_val * dudy;

                // v-momentum convection: u * dv/dx + v * dv/dy
                double dvdx = (v[i+1][j] - v[i-1][j]) / (2*dx);
                double dvdy = (v[i][j+1] - v[i][j-1]) / (2*dy);
                N_v[i][j] = u_val * dvdx + v_val * dvdy;
            }
        }

        if (n == 1) {
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    N_u_old[i][j] = N_u[i][j];
                    N_v_old[i][j] = N_v[i][j];
                }
            }
        }
        //-----------------------------------------------------
        // CN Diffusion
        //-----------------------------------------------------
        cn_diffusion_adi(u, v, N_u, N_u_old, N_v, N_v_old, NX, NY, dx, dy, DT, NU, u_star, v_star);

        for (int j = 0; j < NY; j++) {
            u_star[0][j]     = -u_star[1][j];
            u_star[NX-1][j]  = -u_star[NX-2][j];
            v_star[0][j]     = -v_star[1][j];
            v_star[NX-1][j]  = -v_star[NX-2][j];
        }
        for (int i = 0; i < NX; i++) {
            u_star[i][0]     = -u_star[i][1];
            v_star[i][0]     = -v_star[i][1];
        }
        for (int i = 0; i < NX; i++) {
            u_star[i][NY-1]  = 2.0 * 1.0 - u_star[i][NY-2];
            v_star[i][NY-1]  = -v_star[i][NY-2];
        }
        // =====================================================================


        //-----------------------------------------------------
        // 3) 압력 Poisson 방정식 (Projection)
        //-----------------------------------------------------
        for (int i = 1; i < NX-1; i++) {
            for (int j = 1; j < NY-1; j++) {
                rhs[i][j] = 1.0 / DT 
                            * ( (u_star[i+1][j] - u_star[i-1][j])/(2*dx)
                              + (v_star[i][j+1] - v_star[i][j-1])/(2*dy) );
            }
        }

        for(int i=0; i<NX; i++) for(int j=0; j<NY; j++) p_prime[i][j] = 0.0; // p_prime 초기화
        poisson_pressure_solver(p_prime, rhs, NX, NY, dx, dy, 100000, P_TOL);

        // p_prime 경계조건 적용 (Neumann, dp/dn = 0)
        for (int j = 0; j < NY; j++) {
            p_prime[0][j]     = p_prime[1][j];
            p_prime[NX-1][j]  = p_prime[NX-2][j];
        }
        for (int i = 0; i < NX; i++) {
            p_prime[i][0]     = p_prime[i][1];
            p_prime[i][NY-1]  = p_prime[i][NY-2];
        }

        //-----------------------------------------------------
        // 4) 속도 보정 및 경계조건
        //-----------------------------------------------------
        for (int i = 0; i < NX - 1; i++) {
            for (int j = 1; j < NY - 1; j++) {
                double u_star_face = 0.5 * (u_star[i][j] + u_star[i+1][j]);
                double dp_dx_face = (p_prime[i+1][j] - p_prime[i][j]) / dx;

                u_fc[i][j] = u_star_face 
                        - (DT / RHO) * dp_dx_face;
            }
        }

        for (int i = 1; i < NX - 1; i++) {
            for (int j = 0; j < NY - 1; j++) {
                double v_star_face = 0.5 * (v_star[i][j] + v_star[i][j+1]);
                double dp_dy_face = (p_prime[i][j+1] - p_prime[i][j]) / dy;

                v_fc[i][j] = v_star_face 
                        - (DT / RHO) * dp_dy_face;
            }
        }

        for (int i = 1; i < NX-1; i++) {
            for (int j = 1; j < NY-1; j++) {
                u[i][j] = 0.5 * (u_fc[i][j] + u_fc[i-1][j]);
                v[i][j] = 0.5 * (v_fc[i][j] + v_fc[i][j-1]);
            }
        }

        for (int j = 0; j < NY; j++) {
            u[0][j]     = -u[1][j];
            u[NX-1][j]  = -u[NX-2][j];
            v[0][j]     = -v[1][j]; 
            v[NX-1][j]  = -v[NX-2][j];
        }
        for (int i = 0; i < NX; i++) {
            u[i][0]     = -u[i][1];
            v[i][0]     = -v[i][1];
        }
        // Moving lid (top wall)
        for (int i = 0; i < NX; i++) {
            u[i][NY-1]  = 2.0 * 1.0 - u[i][NY-2];
            v[i][NY-1]  = -v[i][NY-2];
        }

        // 압력 배열 업데이트
        for (int i = 1; i < NX - 1; i++) {
            for (int j = 1; j < NY - 1; j++) {
                p[i][j] += p_prime[i][j];
            }
        }
        for (int j = 0; j < NY; j++) {
            p[0][j]     = p[1][j];
            p[NX-1][j]  = p[NX-2][j];
        }
        for (int i = 0; i < NX; i++) {
            p[i][0]     = p[i][1];
            p[i][NY-1]  = p[i][NY-2];
        }

        //-----------------------------------------------------
        // 5) N_u_old, N_v_old 업데이트 및 수렴 체크
        //-----------------------------------------------------
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                N_u_old[i][j] = N_u[i][j];
                N_v_old[i][j] = N_v[i][j];
            }
        }

        // 수정: 수렴 판정 로직 수정 (정상 상태 도달 시 종료)
        double error = 0.0;
        double all_error = 0.0;
        for (int i = 1; i < NX - 1; i++) {
            for (int j = 1; j < NY - 1; j++) {
                error += fabs(u[i][j] - u_old[i][j]);
                all_error += fabs(u[i][j]);
            }
        }
        error = error / all_error;
        // 발산 감지
        if (isnan(error) || error > 1e5) {
            printf("Divergence detected at step %d with residual %e\n", n, error);
            break;
        }

        print_result(u); // 중간 결과 출력
        printf("Step %d Time %f - residual: %e\n", n, n * DT, error);

        printf("Writing VTK file at step %d...\n", n);
        print_vtk_vector("velocity", n, u, v);
        print_vtk_scalar("pressure", n, p);
    }

    for (int i = 0; i < NX; i++) {
        free(u_star[i]); free(v_star[i]); free(p_prime[i]); free(rhs[i]);
        free(N_u[i]); free(N_u_old[i]); free(N_v[i]); free(N_v_old[i]);
        free(u_old[i]);
    }
    free(u_star); free(v_star); free(p_prime); free(rhs);
    free(N_u); free(N_u_old); free(N_v); free(N_v_old);
    free(u_old);

    for (int i = 0; i < NX - 1; i++) {
        free(u_fc[i]);
    }
    free(u_fc);
    for (int i = 0; i < NX; i++) {
        free(v_fc[i]);
    }
    free(v_fc);
}
