/*
 * This is the main routine of the numerical integrator.
 * It initializes the required memory for its workspace, calls the time-stepping function and controls output.
 */
#include "Verner.h"

void Verner(
    size_t number_of_equations,
    double t,
    double h_step,
    double y[number_of_equations],
    void (*RHS)(double, double *, double *),
    size_t number_of_outputs,
    double t_output[number_of_outputs])
{
    /* Prepare the workspace for the integration. */

    struct workspace *workspace_struct = (struct workspace *)malloc(sizeof(struct workspace));
    workspace_struct->y_temp = (double *)malloc(number_of_equations * sizeof(double));

    for (size_t i = 0; i < NUMBER_OF_STAGES; i++)
        workspace_struct->k[i] = (double *)malloc(number_of_equations * sizeof(double));

    size_t step_counter = 0;                // Counts the steps taken by the integrator so far
    const size_t max_number_of_steps = 500; // Impose a limit on the number of steps to be taken

    while (step_counter < max_number_of_steps)
    {
        /* Make a step. */

        update_stages(number_of_equations, t, h_step, y, RHS, workspace_struct);

        for (size_t i = 0; i < number_of_outputs; i++)
        {
            /* Calculate the auxiliary variable theta such that 0 <= theta <= 1. */
            double theta = (t_output[i] - t) / h_step;

            if (theta >= 0.0 && theta <= 1.0)
                dense_output(number_of_equations, t, theta, h_step, y, workspace_struct);
            else
                continue;
        }

        make_step(number_of_equations, h_step, y, workspace_struct);
        t = t + h_step;
        step_counter = step_counter + 1;
    }
    if (step_counter == max_number_of_steps)
    {
        puts("Reached maximum steps.");
    }

    /* After the integration, clear the workspace. */

    free(workspace_struct->y_temp);
    for (size_t i = 0; i < NUMBER_OF_STAGES; i++)
        free(workspace_struct->k[i]);

    free(workspace_struct);

    return;
}

void update_stages(
    size_t number_of_equations,
    double t,
    double h_step,
    double y[number_of_equations],
    void (*RHS)(double, double *, double *),
    struct workspace *workspace_struct)
{
    /* Allocate space for the stages k1 through k11. Note that k1 is stored in workspace_struct->k[0] due to zero-based indexing. */
    double *const k1 = workspace_struct->k[0];
    double *const k2 = workspace_struct->k[1];
    double *const k3 = workspace_struct->k[2];
    double *const k4 = workspace_struct->k[3];
    double *const k5 = workspace_struct->k[4];
    double *const k6 = workspace_struct->k[5];
    double *const k7 = workspace_struct->k[6];
    double *const k8 = workspace_struct->k[7];
    double *const k9 = workspace_struct->k[8];
    double *const k10 = workspace_struct->k[9];
    double *const k11 = workspace_struct->k[10];

    int i;

    /* k1 stage */
    RHS(t + h_step * RK_c[0], y, k1);

    /* k2 stage */
    for (i = 0; i < number_of_equations; i++)
        workspace_struct->y_temp[i] = y[i] + h_step * (RK_A[1][0] * k1[i]);
    RHS(t + h_step * RK_c[1], workspace_struct->y_temp, k2);

    /* k3 stage */
    for (i = 0; i < number_of_equations; i++)
        workspace_struct->y_temp[i] = y[i] + h_step * (RK_A[2][0] * k1[i] + RK_A[2][1] * k2[i]);
    RHS(t + h_step * RK_c[2], workspace_struct->y_temp, k3);

    /* k4 stage */
    for (i = 0; i < number_of_equations; i++)
        workspace_struct->y_temp[i] = y[i] + h_step * (RK_A[3][0] * k1[i] + RK_A[3][1] * k2[i] + RK_A[3][2] * k3[i]);
    RHS(t + h_step * RK_c[3], workspace_struct->y_temp, k4);

    /* k5 stage */
    for (i = 0; i < number_of_equations; i++)
        workspace_struct->y_temp[i] = y[i] + h_step * (RK_A[4][0] * k1[i] + RK_A[4][1] * k2[i] + RK_A[4][2] * k3[i] + RK_A[4][3] * k4[i]);
    RHS(t + h_step * RK_c[4], workspace_struct->y_temp, k5);

    /* k6 stage */
    for (i = 0; i < number_of_equations; i++)
        workspace_struct->y_temp[i] = y[i] + h_step * (RK_A[5][0] * k1[i] + RK_A[5][1] * k2[i] + RK_A[5][2] * k3[i] + RK_A[5][3] * k4[i] + RK_A[5][4] * k5[i]);
    RHS(t + h_step * RK_c[5], workspace_struct->y_temp, k6);

    /* k7 stage */
    for (i = 0; i < number_of_equations; i++)
        workspace_struct->y_temp[i] = y[i] + h_step * (RK_A[6][0] * k1[i] + RK_A[6][1] * k2[i] + RK_A[6][2] * k3[i] + RK_A[6][3] * k4[i] + RK_A[6][4] * k5[i] + RK_A[6][5] * k6[i]);
    RHS(t + h_step * RK_c[6], workspace_struct->y_temp, k7);

    /* k8 stage */
    for (i = 0; i < number_of_equations; i++)
        workspace_struct->y_temp[i] = y[i] + h_step * (RK_A[7][0] * k1[i] + RK_A[7][1] * k2[i] + RK_A[7][2] * k3[i] + RK_A[7][3] * k4[i] + RK_A[7][4] * k5[i] + RK_A[7][5] * k6[i] + RK_A[7][6] * k7[i]);
    RHS(t + h_step * RK_c[7], workspace_struct->y_temp, k8);

    /* k9 stage */
    for (i = 0; i < number_of_equations; i++)
        workspace_struct->y_temp[i] = y[i] + h_step * (RK_A[8][0] * k1[i] + RK_A[8][1] * k2[i] + RK_A[8][2] * k3[i] + RK_A[8][3] * k4[i] + RK_A[8][4] * k5[i] + RK_A[8][5] * k6[i] + RK_A[8][6] * k7[i] + RK_A[8][7] * k8[i]);
    RHS(t + h_step * RK_c[5], workspace_struct->y_temp, k9);

    /* k10 stage */
    for (i = 0; i < number_of_equations; i++)
        workspace_struct->y_temp[i] = y[i] + h_step * (RK_A[9][0] * k1[i] + RK_A[9][1] * k2[i] + RK_A[9][2] * k3[i] + RK_A[9][3] * k4[i] + RK_A[9][4] * k5[i] + RK_A[9][5] * k6[i] + RK_A[9][6] * k7[i] + RK_A[9][7] * k8[i] + RK_A[9][8] * k9[i]);
    RHS(t + h_step * RK_c[6], workspace_struct->y_temp, k10);

    /* k11 stage */
    for (i = 0; i < number_of_equations; i++)
        workspace_struct->y_temp[i] = y[i] + h_step * (RK_A[10][0] * k1[i] + RK_A[10][1] * k2[i] + RK_A[10][2] * k3[i] + RK_A[10][3] * k4[i] + RK_A[10][4] * k5[i] + RK_A[10][5] * k6[i] + RK_A[10][6] * k7[i] + RK_A[10][7] * k8[i] + RK_A[10][8] * k9[i] + RK_A[10][9] * k10[i]);
    RHS(t + h_step * RK_c[7], workspace_struct->y_temp, k11);

    return;
}

void make_step(int number_of_equations, double h_step, double y[number_of_equations], struct workspace *workspace_struct)
{
    /* Fetch the stages k1 and k8 through k11. */
    double *const k1 = workspace_struct->k[0];
    double *const k8 = workspace_struct->k[7];
    double *const k9 = workspace_struct->k[8];
    double *const k10 = workspace_struct->k[9];
    double *const k11 = workspace_struct->k[10];

    /* Combine the k-stages in the final sum to forward the solution. */
    for (int i = 0; i < number_of_equations; i++)
        y[i] = y[i] + h_step * (RK_b[0] * k1[i] + RK_b[7] * k8[i] + RK_b[8] * k9[i] + RK_b[9] * k10[i] + RK_b[10] * k11[i]);

    return;
}
