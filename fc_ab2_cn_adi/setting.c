#include "setting.h"

void print_result(double **solution) {
    int n_stride = 8;
    int x_stride = (NX > n_stride) ? NX / n_stride : 1;
    int y_stride = (NY > n_stride) ? NY / n_stride : 1;
    for (int j = NY - 1; j >= 0; j -= y_stride) {
        for (int i = 0; i < NX; i += x_stride) {
            printf("%.6f\t", solution[i][j]);
        }
        printf("\n");
    }
}

void print_file(double **solution, int number) {
    char file_name[256];
    sprintf(file_name, "solution_output_%d.txt", number);
    FILE *fp = fopen(file_name, "w");
    if (fp == NULL) {
        perror("Error opening output file");
        return;
    }
    fprintf(fp, "%d %d\n", NX, NY);

    for (int j = 0; j < NY; j++) {
        for (int i = 0; i < NX; i++) {
            fprintf(fp, "%lf ", solution[j][i]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}

void print_vtk_file(double **field, int step, const char *varname) {
    char file_name[256];
    sprintf(file_name, "outputs/%s_t%04d.vtk", varname, step);
    FILE *fp = fopen(file_name, "w");
    if (fp == NULL) {
        perror("Error opening output file");
        return;
    }
    
    const int pNX = NX - 2;
    const int pNY = NY - 2;
    
    // 물리적 격자 간격
    const double dx = 1.0 / (pNX - 1);
    const double dy = 1.0 / (pNY - 1);

    // VTK header
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "%s field (physical domain), step %d\n", varname, step);
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET STRUCTURED_POINTS\n");
    
    // 물리적 격자점 수로 DIMENSIONS 설정
    fprintf(fp, "DIMENSIONS %d %d 1\n", pNX, pNY);
    
    // 물리적 원점과 격자 간격 설정
    fprintf(fp, "ORIGIN 0.0 0.0 0.0\n");
    fprintf(fp, "SPACING %e %e 1.0\n", dx, dy);
    
    // 데이터 속성 (물리적 점의 개수)
    fprintf(fp, "POINT_DATA %d\n", pNX * pNY);
    fprintf(fp, "SCALARS %s double 1\n", varname);
    fprintf(fp, "LOOKUP_TABLE default\n");

    // 데이터 출력 (물리적 도메인만)
    for (int j = 1; j < NY - 1; j++) {
        for (int i = 1; i < NX - 1; i++) {
            fprintf(fp, "%e\n", field[i][j]);
        }
    }
    
    fclose(fp);
}

void print_vtk_scalar(const char* varname, int step, double** field) {
    char file_name[256];
    // outputs 디렉토리 미리 생성
    sprintf(file_name, "outputs/%s_Re%d_t%04d.vtk", varname, (int)_CONST_Re, step);
    FILE *fp = fopen(file_name, "w");
    if (fp == NULL) {
        perror("Error: Cannot open VTK scalar file for writing. Ensure 'outputs' directory exists.");
        return;
    }

    const int pNX = NX - 2;
    const int pNY = NY - 2;
    
    const double dx = 1.0 / (pNX - 1);
    const double dy = 1.0 / (pNY - 1);

    // --- VTK Header ---
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "%s field, step %d\n", varname, step);
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET STRUCTURED_POINTS\n");
    
    // 물리적 격자점 수로 DIMENSIONS 설정
    fprintf(fp, "DIMENSIONS %d %d 1\n", pNX, pNY);
    
    // 물리적 원점과 격자 간격 설정
    fprintf(fp, "ORIGIN 0.0 0.0 0.0\n");
    fprintf(fp, "SPACING %e %e 1.0\n", dx, dy);
    
    // 데이터 속성 (물리적 점의 개수)
    fprintf(fp, "POINT_DATA %d\n", pNX * pNY);
    fprintf(fp, "SCALARS %s double 1\n", varname);
    fprintf(fp, "LOOKUP_TABLE default\n");

    // 데이터 출력 (물리적 도메인만, Ghost cell 제외)
    for (int j = 1; j < NY - 1; j++) {
        for (int i = 1; i < NX - 1; i++) {
            fprintf(fp, "%e\n", field[i][j]);
        }
    }
    
    fclose(fp);
}

void print_vtk_vector(const char* varname, int step, double** u, double** v) {
    char file_name[256];
    sprintf(file_name, "outputs/%s_Re%d_t%04d.vtk", varname, (int)_CONST_Re, step);
    FILE *fp = fopen(file_name, "w");
    if (fp == NULL) {
        perror("Error: Cannot open VTK vector file for writing. Ensure 'outputs' directory exists.");
        return;
    }

    const int pNX = NX - 2;
    const int pNY = NY - 2;
    
    const double dx = 1.0 / (pNX - 1);
    const double dy = 1.0 / (pNY - 1);

    // --- VTK Header ---
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "%s field, step %d\n", varname, step);
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET STRUCTURED_POINTS\n");
    
    fprintf(fp, "DIMENSIONS %d %d 1\n", pNX, pNY);
    fprintf(fp, "ORIGIN 0.0 0.0 0.0\n");
    fprintf(fp, "SPACING %e %e 1.0\n", dx, dy);
    
    fprintf(fp, "POINT_DATA %d\n", pNX * pNY);

    fprintf(fp, "VECTORS %s double\n", varname);

    for (int j = 1; j < NY - 1; j++) {
        for (int i = 1; i < NX - 1; i++) {
            // z-성분은 0.0으로 설정
            fprintf(fp, "%e %e 0.0\n", u[i][j], v[i][j]);
        }
    }
    
    fclose(fp);
}

void print_both(FILE *fp, const char *format, ...) {
    va_list args;
    va_start(args, format);
    vprintf(format, args);
    va_end(args);
    
    va_start(args, format);
    vfprintf(fp, format, args);
    va_end(args);
}

void print_current_time(FILE *fp, const char *prefix) {
    struct timeval tv;
    gettimeofday(&tv, NULL);

    struct tm *tm_info = localtime(&tv.tv_sec);
    char buffer[64];
    strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", tm_info);

    // 밀리초 붙이기
    // fprintf(fp, "%s%s.%03ld\n", prefix, buffer, tv.tv_usec / 1000);
    // printf("%s%s.%03ld\n", prefix, buffer, tv.tv_usec / 1000);
}


