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
void Verner(size_t number_of_equations, double t, double h_step, double y[number_of_equations], void (*RHS)(double, double *, double *), size_t number_of_outputs, double t_output[number_of_outputs]);

void update_stages(size_t number_of_equations, double t, double h_step, double y[number_of_equations], void (*RHS)(double, double *, double *), struct workspace *workspace_struct);
void make_step(int number_of_equations, double h_step, double y[number_of_equations], struct workspace *workspace_struct);

// dense_output.c
void dense_output(size_t number_of_equations, double t, double theta, double h_step, double y[number_of_equations], struct workspace *workspace_struct);
