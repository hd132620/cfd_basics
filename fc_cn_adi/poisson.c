#include "poisson.h"

static void tdma(int n, double *a, double *b, double *c, double *d, double *x) {
    double *cp = (double*)malloc(n * sizeof(double));
    double *dp = (double*)malloc(n * sizeof(double));
    
    // Forward elimination
    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];
    for (int i = 1; i < n; i++) {
        double m = b[i] - a[i] * cp[i-1];
        cp[i] = c[i] / m;
        dp[i] = (d[i] - a[i] * dp[i-1]) / m;
    }
    
    // Backward substitution
    x[n-1] = dp[n-1];
    for (int i = n-2; i >= 0; i--) {
        x[i] = dp[i] - cp[i] * x[i+1];
    }
    
    free(cp); free(dp);
}

void poisson_pressure_solver(
    double **p,
    double **rhs,
    int     nx,
    int     ny,
    double  dx,
    double  dy,
    int     max_iter,
    double  tol
) {
    double alpha = dx*dx/(dy*dy);
    double beta  = dy*dy/(dx*dx);
    int   maxn   = (nx>ny? nx: ny);

    // tridiagonal buffers
    double *a   = calloc(maxn,sizeof(double));
    double *b   = calloc(maxn,sizeof(double));
    double *c   = calloc(maxn,sizeof(double));
    double *d   = calloc(maxn,sizeof(double));
    double *sol = calloc(maxn,sizeof(double));

    // Initial Residual
    double initial_max_res = 0.0;
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            double lapx = (p[i+1][j] - 2*p[i][j] + p[i-1][j])/(dx*dx);
            double lapy = (p[i][j+1] - 2*p[i][j] + p[i][j-1])/(dy*dy);
            double res  = fabs(lapx + lapy - rhs[i][j]);
            if (res > initial_max_res) initial_max_res = res;
        }
    }

    for (int it = 0; it < max_iter; it++) {
        // x-sweep
        for (int j = 1; j < ny-1; j++) {
            for (int i = 0; i < nx; i++) {
                if (i == 0) {
                    b[0] =  1.0; c[0] = -1.0; d[0] =  0.0;
                    a[0] =  0.0;
                } else if (i == nx-1) {
                    a[nx-1] = -1.0;
                    b[nx-1] =  1.0;
                    c[nx-1] =  0.0;
                    d[nx-1] =  0.0;
                } else {
                    a[i] = 1.0;
                    b[i] = -2.0*(1.0 + alpha);
                    c[i] = 1.0;
                    d[i] = rhs[i][j]*dx*dx
                         - alpha*(p[i][j+1] + p[i][j-1]);
                }
            }
            tdma(nx, a, b, c, d, sol);
            for (int i = 0; i < nx; i++) {
                p[i][j] = sol[i];
            }
        }

        // y-sweep
        for (int i = 1; i < nx-1; i++) {
            for (int j = 0; j < ny; j++) {
                if (j == 0) {
                    b[0] =  1.0; c[0] = -1.0; d[0] =  0.0;
                    a[0] =  0.0; // 사용되지 않음
                } else if (j == ny-1) {
                    a[ny-1] = -1.0;
                    b[ny-1] =  1.0;
                    c[ny-1] =  0.0;
                    d[ny-1] =  0.0;
                } else {
                    a[j] = 1.0;
                    b[j] = -2.0*(1.0 + beta);
                    c[j] = 1.0;
                    d[j] = rhs[i][j]*dy*dy
                         - beta*(p[i+1][j] + p[i-1][j]);
                }
            }
            tdma(ny, a, b, c, d, sol);
            for (int j = 0; j < ny; j++) {
                p[i][j] = sol[j];
            }
        }

        // Current Residual
        double max_res = 0.0;
        for (int i = 1; i < nx-1; i++) {
            for (int j = 1; j < ny-1; j++) {
                double lapx = (p[i+1][j] - 2*p[i][j] + p[i-1][j])/(dx*dx);
                double lapy = (p[i][j+1] - 2*p[i][j] + p[i][j-1])/(dy*dy);
                double res  = fabs(lapx + lapy - rhs[i][j]);
                if (res > max_res) max_res = res;
            }
        }

        // printf("Poisson Iteration %d: max residual = %e\n", it, max_res);
        
        // 상대 잔차로 수렴 판정
        if (max_res / initial_max_res < tol) {
            printf("Poisson converged at iteration %d with relative residual %e\n", it, max_res / initial_max_res);
            break;
        }

        // 만약 최대 반복 횟수에 도달했다면 경고 메시지 출력
        if (it == max_iter - 1) {
            printf("WARNING: Poisson solver did not converge. Final relative residual: %e\n", max_res / initial_max_res);
        }
    }

    free(a); free(b); free(c); free(d); free(sol);
}


static void cn_x_sweep(
    double **phi_old, double **N_phi, double **N_phi_old,
    double dx, double dy, double dt, double nu, int nx, int ny,
    double **phi_tilde
) {
    int M = nx - 2;
    if (M <= 0) return;

    double *a = (double*)malloc(M * sizeof(double));
    double *b = (double*)malloc(M * sizeof(double));
    double *c = (double*)malloc(M * sizeof(double));
    double *d = (double*)malloc(M * sizeof(double));
    double *sol = (double*)malloc(M * sizeof(double));

    double alpha_x = nu * dt / (2.0 * dx * dx);
    double alpha_y = nu * dt / (2.0 * dy * dy);

    for (int j = 1; j < ny - 1; j++) {
        for (int i = 0; i < M; i++) {
            a[i] = -alpha_x;
            b[i] = 1.0 + 2.0 * alpha_x;
            c[i] = -alpha_x;
        }

        // Calculate RHS of the equation
        for (int i_idx = 0; i_idx < M; i_idx++) {
            int i = i_idx + 1;
            double conv_ab2 = 0.5 * (3.0 * N_phi[i][j] - N_phi_old[i][j]);
            double lap_y_old = (phi_old[i][j+1] - 2*phi_old[i][j] + phi_old[i][j-1]) / (dy*dy);

            d[i_idx] = phi_old[i][j] - dt * conv_ab2 + 0.5 * dt * nu * lap_y_old;
        }

        tdma(M, a, b, c, d, sol);

        for (int i_idx = 0; i_idx < M; i_idx++) {
            phi_tilde[i_idx + 1][j] = sol[i_idx];
        }
    }
    free(a); free(b); free(c); free(d); free(sol);
}

static void cn_y_sweep(
    double **phi_old, double **phi_tilde,
    double dx, double dy, double dt, double nu, int nx, int ny,
    double **phi_new
) {
    int M = ny - 2;
    if (M <= 0) return;

    double *a = (double*)malloc(M * sizeof(double));
    double *b = (double*)malloc(M * sizeof(double));
    double *c = (double*)malloc(M * sizeof(double));
    double *d = (double*)malloc(M * sizeof(double));
    double *sol = (double*)malloc(M * sizeof(double));

    double alpha_x = nu * dt / (2.0 * dx * dx);
    double alpha_y = nu * dt / (2.0 * dy * dy);
    
    for (int i = 1; i < nx - 1; i++) {
        for (int j = 0; j < M; j++) {
            a[j] = -alpha_y;
            b[j] = 1.0 + 2.0 * alpha_y;
            c[j] = -alpha_y;
        }

        for (int j_idx = 0; j_idx < M; j_idx++) {
            int j = j_idx + 1;
            double lap_x_tilde = (phi_tilde[i+1][j+1] - 2.0*phi_tilde[i][j] + phi_tilde[i-1][j]) / (dx*dx);
            d[j_idx] = phi_tilde[i][j] + 0.5 * dt * nu * lap_x_tilde;
        }

        tdma(M, a, b, c, d, sol);

        for (int j_idx = 0; j_idx < M; j_idx++) {
            phi_new[i][j_idx + 1] = sol[j_idx];
        }
    }
    free(a); free(b); free(c); free(d); free(sol);
}



// --- Main ADI Function ---
void cn_diffusion_adi(
    double **u_old, double **v_old,
    double **N_u,   double **N_u_old, 
    double **N_v,   double **N_v_old,
    int nx, int ny, double dx, double dy,
    double dt, double nu,
    double **u_star, double **v_star // Output
) {
    // Allocate temporary fields for intermediate results
    double **u_tilde = (double**)malloc(nx * sizeof(double*));
    double **v_tilde = (double**)malloc(nx * sizeof(double*));
    for (int i = 0; i < nx; i++) {
        u_tilde[i] = (double*)calloc(ny, sizeof(double));
        v_tilde[i] = (double*)calloc(ny, sizeof(double));
    }

    // --- Process u-momentum ---
    // 1. Predictor step (implicit in x)
    cn_x_sweep(u_old, N_u, N_u_old, dx, dy, dt, nu, nx, ny, u_tilde);
    // 2. Corrector step (implicit in y)
    cn_y_sweep(u_old, u_tilde, dx, dy, dt, nu, nx, ny, u_star);

    // --- Process v-momentum ---
    // 1. Predictor step (implicit in x)
    cn_x_sweep(v_old, N_v, N_v_old, dx, dy, dt, nu, nx, ny, v_tilde);
    // 2. Corrector step (implicit in y)
    cn_y_sweep(v_old, v_tilde, dx, dy, dt, nu, nx, ny, v_star);

    for (int i = 0; i < nx; i++) {
        free(u_tilde[i]);
        free(v_tilde[i]);
    }
    free(u_tilde);
    free(v_tilde);
}