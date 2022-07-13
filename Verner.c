/*
 * This is the main routine of the numerical integrator.
 * It initializes the required memory for its workspace, calls the time-stepping function and controls output.
 */
#include "Verner.h"

void Verner(void (*RHS)(double, double *, double *), int n_equat, double t, double y[n_equat], double h_step, double t_output[])
{
    /* Prepare the workspace for the integration. */

    struct workspace *myWorkspace = (struct workspace *)malloc(sizeof(struct workspace));
    myWorkspace->y_temp = (double *)malloc(n_equat * sizeof(double));

    for (size_t i = 0; i < NUMBER_OF_STAGES; i++)
        myWorkspace->k[i] = (double *)malloc(n_equat * sizeof(double));

    int stepsCounter = 0; // A variable that counts the steps taken by the integrator so far

    while (stepsCounter < 500)
    {
        /* Make a step. */

        UpdateStages(RHS, n_equat, t, y, h_step, myWorkspace);
        DenseOutput(n_equat, t, y, h_step, t_output, myWorkspace);
        Step(n_equat, t, y, h_step, myWorkspace);
        t = t + h_step;
        stepsCounter = stepsCounter + 1;
    }
    if (stepsCounter == 500)
    {
        puts("Reached maximum steps.");
    }

    /* After the integration, clear the workspace. */

    free(myWorkspace->y_temp);
    for (size_t i = 0; i < NUMBER_OF_STAGES; i++)
        free(myWorkspace->k[i]);

    free(myWorkspace);

    return;
}

/* This module forms the integrator core that advances the numerical solution. */
void UpdateStages(void (*RHS)(double, double *, double *), int n_equat, double t, double y[n_equat], double h_step, struct workspace *myWorkspace)
{
    /* Allocate space for the stages k1 through k11. Note that k1 is stored in myWorkspace->k[0] due to zero-based indexing. */
    double *const k1 = myWorkspace->k[0];
    double *const k2 = myWorkspace->k[1];
    double *const k3 = myWorkspace->k[2];
    double *const k4 = myWorkspace->k[3];
    double *const k5 = myWorkspace->k[4];
    double *const k6 = myWorkspace->k[5];
    double *const k7 = myWorkspace->k[6];
    double *const k8 = myWorkspace->k[7];
    double *const k9 = myWorkspace->k[8];
    double *const k10 = myWorkspace->k[9];
    double *const k11 = myWorkspace->k[10];

    int i;

    /* k1 stage */
    RHS(t + h_step * RK_c[0], y, k1);

    /* k2 stage */
    for (i = 0; i < n_equat; i++)
        myWorkspace->y_temp[i] = y[i] + h_step * (RK_A[1][0] * k1[i]);
    RHS(t + h_step * RK_c[1], myWorkspace->y_temp, k2);

    /* k3 stage */
    for (i = 0; i < n_equat; i++)
        myWorkspace->y_temp[i] = y[i] + h_step * (RK_A[2][0] * k1[i] + RK_A[2][1] * k2[i]);
    RHS(t + h_step * RK_c[2], myWorkspace->y_temp, k3);

    /* k4 stage */
    for (i = 0; i < n_equat; i++)
        myWorkspace->y_temp[i] = y[i] + h_step * (RK_A[3][0] * k1[i] + RK_A[3][1] * k2[i] + RK_A[3][2] * k3[i]);
    RHS(t + h_step * RK_c[3], myWorkspace->y_temp, k4);

    /* k5 stage */
    for (i = 0; i < n_equat; i++)
        myWorkspace->y_temp[i] = y[i] + h_step * (RK_A[4][0] * k1[i] + RK_A[4][1] * k2[i] + RK_A[4][2] * k3[i] + RK_A[4][3] * k4[i]);
    RHS(t + h_step * RK_c[4], myWorkspace->y_temp, k5);

    /* k6 stage */
    for (i = 0; i < n_equat; i++)
        myWorkspace->y_temp[i] = y[i] + h_step * (RK_A[5][0] * k1[i] + RK_A[5][1] * k2[i] + RK_A[5][2] * k3[i] + RK_A[5][3] * k4[i] + RK_A[5][4] * k5[i]);
    RHS(t + h_step * RK_c[5], myWorkspace->y_temp, k6);

    /* k7 stage */
    for (i = 0; i < n_equat; i++)
        myWorkspace->y_temp[i] = y[i] + h_step * (RK_A[6][0] * k1[i] + RK_A[6][1] * k2[i] + RK_A[6][2] * k3[i] + RK_A[6][3] * k4[i] + RK_A[6][4] * k5[i] + RK_A[6][5] * k6[i]);
    RHS(t + h_step * RK_c[6], myWorkspace->y_temp, k7);

    /* k8 stage */
    for (i = 0; i < n_equat; i++)
        myWorkspace->y_temp[i] = y[i] + h_step * (RK_A[7][0] * k1[i] + RK_A[7][1] * k2[i] + RK_A[7][2] * k3[i] + RK_A[7][3] * k4[i] + RK_A[7][4] * k5[i] + RK_A[7][5] * k6[i] + RK_A[7][6] * k7[i]);
    RHS(t + h_step * RK_c[7], myWorkspace->y_temp, k8);

    /* k9 stage */
    for (i = 0; i < n_equat; i++)
        myWorkspace->y_temp[i] = y[i] + h_step * (RK_A[8][0] * k1[i] + RK_A[8][1] * k2[i] + RK_A[8][2] * k3[i] + RK_A[8][3] * k4[i] + RK_A[8][4] * k5[i] + RK_A[8][5] * k6[i] + RK_A[8][6] * k7[i] + RK_A[8][7] * k8[i]);
    RHS(t + h_step * RK_c[5], myWorkspace->y_temp, k9);

    /* k10 stage */
    for (i = 0; i < n_equat; i++)
        myWorkspace->y_temp[i] = y[i] + h_step * (RK_A[9][0] * k1[i] + RK_A[9][1] * k2[i] + RK_A[9][2] * k3[i] + RK_A[9][3] * k4[i] + RK_A[9][4] * k5[i] + RK_A[9][5] * k6[i] + RK_A[9][6] * k7[i] + RK_A[9][7] * k8[i] + RK_A[9][8] * k9[i]);
    RHS(t + h_step * RK_c[6], myWorkspace->y_temp, k10);

    /* k11 stage */
    for (i = 0; i < n_equat; i++)
        myWorkspace->y_temp[i] = y[i] + h_step * (RK_A[10][0] * k1[i] + RK_A[10][1] * k2[i] + RK_A[10][2] * k3[i] + RK_A[10][3] * k4[i] + RK_A[10][4] * k5[i] + RK_A[10][5] * k6[i] + RK_A[10][6] * k7[i] + RK_A[10][7] * k8[i] + RK_A[10][8] * k9[i] + RK_A[10][9] * k10[i]);
    RHS(t + h_step * RK_c[7], myWorkspace->y_temp, k11);

    return;
}

void Step(int n_equat, double t, double y[n_equat], double h_step, struct workspace *myWorkspace)
{
    /* Fetch the stages k1 and k8 through k11. */
    double *const k1 = myWorkspace->k[0];
    double *const k8 = myWorkspace->k[7];
    double *const k9 = myWorkspace->k[8];
    double *const k10 = myWorkspace->k[9];
    double *const k11 = myWorkspace->k[10];

    /* Combine the k-stages in the final sum to forward the solution. */
    for (int i = 0; i < n_equat; i++)
        y[i] = y[i] + h_step * (RK_b[0] * k1[i] + RK_b[7] * k8[i] + RK_b[8] * k9[i] + RK_b[9] * k10[i] + RK_b[10] * k11[i]);

    return;
}
