/* Cooper-Verner 11 stage order 8 Runge-Kutta scheme */

#include <stdio.h>
#include <stdlib.h>

/* Preprocessor global variable declarations */
#define ORDER 8
#define NUMBER_OF_STAGES 11
#define sqrt21 4.582575694955840006588047193728

/* Coefficients arrays for the scheme */
extern const double RK_c[NUMBER_OF_STAGES], RK_b[NUMBER_OF_STAGES];
extern const double RK_A[NUMBER_OF_STAGES][NUMBER_OF_STAGES];

/* Workspace struct declaration */
struct workspace
{
    double *k[NUMBER_OF_STAGES];
    double *y_temp;
};

/* Function prototypes */
// Verner.c
void Verner(void (*RHS)(double, double *, double *), int n_equat, double t, double y[n_equat], double h_step, double t_output[]);
void UpdateStages(void (*RHS)(double, double *, double *), int n_equat, double t, double y[n_equat], double h_step, struct workspace *myWorkspace);
void DenseOutput(int n_equat, double t, double y[n_equat], double h_step, double t_output[], struct workspace *myWorkspace);
void Step(int n_equat, double t, double y[n_equat], double h_step, struct workspace *myWorkspace);
