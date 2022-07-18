/*
 * This is the main function that calls the IMEX numerical integrator routine for execution.
 * The function to be integrated and its initial conditions, as well as the step-size, are provided by the user.
 */
#include "tests.h"
#include "Verner.h"

int main(void)
{
    puts("Program: Cooper-Verner 11-stage 8th-order Runge-Kutta Scheme");
    printf("\n");

    /* Specify the ODE function or system to be solved. */

    void (*RHS)(double, double *, double *);
    RHS = RHSCosine;
    size_t number_of_equations = 1;

    /* Allocate memory to the global integration variables. */

    double t;

    double *y;
    y = (double *)malloc(number_of_equations * sizeof(double));

    double *dydt;
    dydt = (double *)malloc(number_of_equations * sizeof(double));

    /* Set the initial conditions of the ODE system. */

    t = 0.0;
    y[0] = 0.0;

    double h_step = 0.1; // The constant step size taken by the numerical integrator
    printf("h_step = %G\n", h_step);
    printf("\n");

    /* Specify the points at which to query the integrator. */

    size_t number_of_outputs = 16;
    double t_output[number_of_outputs];
    for (int i = 0; i <= 15; i++)
        t_output[i] = (double)i;

    /* Issue a call to the primary integrator method. */

    Verner(number_of_equations, t, h_step, y, RHS, number_of_outputs, t_output);

    /* Once the integrator is done, free memory of the global variables. */

    free(y);
    free(dydt);

    printf("Press ENTER key to end.\n");
    getchar();
    return 0;
}
