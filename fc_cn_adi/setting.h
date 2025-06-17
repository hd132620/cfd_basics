#ifndef __SETTING_H__
#define __SETTING_H__

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>
#include <float.h>

#include "navier_stokes.h"

#define DT 5e-3
#define NT 10000
#define NX 64
#define NY 64
#define eps 1e-6
#define P_TOL 1e-6

#define NXNY (NX*NY)
#define ONE_AXIS(X,Y) (((Y)*(NX))+(X))

double ** solution;

void print_result(double ** solution);
void print_file(double ** solution, int number);
void print_both(FILE *fp, const char *format, ...);
void print_current_time(FILE *fp, const char *prefix);
void print_vtk_file(double **field, int number, const char *varname);
void print_vtk_scalar(const char* varname, int step, double** field);
void print_vtk_vector(const char* varname, int step, double** u, double** v);

#endif
