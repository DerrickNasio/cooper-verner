/*
 * main.c
 * Demonstration program
 */
#include <stdlib.h>

#include "tests.h"
#include "verner.h"

int main(void)
{
    puts("Program: Cooper-Verner 11-stage 8th-order Runge-Kutta Scheme");
    printf("\n");

    /* Specify the ODE function or system to be solved. */

    void (*RHS)(double, double *, double *);
    RHS = RHSSineCosine;
    size_t number_of_equations = 2;

    /* Allocate memory to the global integration variables. */

    double t;

    double *y;
    y = (double *)malloc(number_of_equations * sizeof(double));

    double *dydt;
    dydt = (double *)malloc(number_of_equations * sizeof(double));

    /* Set the initial conditions of the ODE system. */

    t = 0.0;
    y[0] = 1.0;
    y[1] = 0.0;

    double h_step = 0.1; // The constant step size taken by the numerical integrator
    printf("h_step = %G\n", h_step);
    printf("\n");

    /* Specify the points at which to query the integrator. */
    size_t number_of_outputs = 13;
    double t_output[number_of_outputs];
    for (size_t i = 0; i < number_of_outputs; i++)
        t_output[i] = ((double)i / 12.0) * 3.14159265358979323846;

    /* Issue a call to the primary integrator method. */
    verner(number_of_equations, t, h_step, y, RHS, number_of_outputs, t_output);

    /* Print the local time. */
    printf("\n");
    timestamp();

    /* Once the integrator is done, free memory of the global variables. */
    free(y);
    free(dydt);

    printf("Press ENTER key to end.\n");
    getchar();
    return 0;
}

void verner(
    size_t number_of_equations,
    double t,
    double h_step,
    double y[number_of_equations],
    void (*RHS)(double, double *, double *),
    size_t number_of_outputs,
    double t_output[number_of_outputs])
{
    /* Prepare the workspace for the integration. */

    verner_state_t *state;
    state = verner_alloc(number_of_equations);

    size_t step_counter = 0;                // Counts the steps taken by the integrator so far
    const size_t max_number_of_steps = 500; // Impose a limit on the number of steps to be taken

    while (step_counter < max_number_of_steps)
    {
        /* Make a step. */

        verner_compute_stages(state, number_of_equations, t, h_step, y, RHS);

        for (size_t i = 0; i < number_of_outputs; i++)
        {
            /* Calculate the auxiliary variable theta such that 0 <= theta <= 1. */
            double theta = (t_output[i] - t) / h_step;

            if (theta >= 0.0 && theta <= 1.0)
                verner_dense_output(state, number_of_equations, t, theta, h_step, y);
            else
                continue;
        }

        verner_apply(state, number_of_equations, h_step, y);
        t = t + h_step;
        step_counter = step_counter + 1;
    }
    if (step_counter == max_number_of_steps)
    {
        puts("Reached maximum steps.");
    }

    /* After the integration, clear the workspace. */

    verner_free(state);

    return;
}
