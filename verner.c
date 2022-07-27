/*
 * verner.c
 * Runge-Kutta, 8th-order, Cooper-Verner
 *
 * Some Explicit Runge-Kutta Methods of High Order
 * Cooper, G. J., and Verner, J. H.
 * SIAM Journal on Numerical Analysis, vol. 9, no. 3, September 1972, pp. 389–405.
 */
#include <stdio.h>
#include <stdlib.h>

#include "verner.h"

// Nodes
static const double RK_c[NUMBER_OF_STAGES] = {
    0.0,
    1.0 / 2.0,
    1.0 / 2.0,
    (7.0 + sqrt21) / 14.0,
    (7.0 + sqrt21) / 14.0,
    1.0 / 2.0,
    (7.0 - sqrt21) / 14.0,
    (7.0 - sqrt21) / 14.0,
    1.0 / 2.0,
    (7.0 + sqrt21) / 14.0,
    1.0};

// Stage weights for final sum
static const double RK_b[NUMBER_OF_STAGES] = {
    1.0 / 20.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    49.0 / 180.0,
    16.0 / 45.0,
    49.0 / 180.0,
    1.0 / 20.0};

// Coupling coefficients
static const double RK_A[NUMBER_OF_STAGES][NUMBER_OF_STAGES] = {
    {0.0},
    {1.0 / 2.0},
    {1.0 / 4.0,
     1.0 / 4.0},
    {1.0 / 7.0,
     -(7.0 + 3.0 * sqrt21) / 98.0,
     (21.0 + 5.0 * sqrt21) / 49.0},
    {(11.0 + sqrt21) / 84.0,
     0.0,
     (18.0 + 4.0 * sqrt21) / 63.0,
     (21.0 - sqrt21) / 252.0},
    {(5.0 + sqrt21) / 48.0,
     0.0,
     (9.0 + sqrt21) / 36.0,
     (-231.0 + 14.0 * sqrt21) / 360.0,
     (63.0 - 7.0 * sqrt21) / 80.0},
    {(10.0 - sqrt21) / 42.0,
     0.0,
     (-432.0 + 92.0 * sqrt21) / 315.0,
     (633.0 - 145.0 * sqrt21) / 90.0,
     (-504.0 + 115.0 * sqrt21) / 70.0,
     (63.0 - 13.0 * sqrt21) / 35.0},
    {1.0 / 14.0,
     0.0,
     0.0,
     0.0,
     (14.0 - 3.0 * sqrt21) / 126.0,
     (13.0 - 3.0 * sqrt21) / 63.0,
     1.0 / 9.0},
    {1.0 / 32.0,
     0.0,
     0.0,
     0.0,
     (91.0 - 21.0 * sqrt21) / 576.0,
     11.0 / 72.0,
     -(385.0 + 75.0 * sqrt21) / 1152.0,
     (63.0 + 13.0 * sqrt21) / 128.0},
    {1.0 / 14.0,
     0.0,
     0.0,
     0.0,
     1.0 / 9.0,
     -(733.0 + 147.0 * sqrt21) / 2205.0,
     (515.0 + 111.0 * sqrt21) / 504.0,
     -(51.0 + 11.0 * sqrt21) / 56.0,
     (132.0 + 28.0 * sqrt21) / 245.0},
    {0.0,
     0.0,
     0.0,
     0.0,
     (-42.0 + 7.0 * sqrt21) / 18.0,
     (-18.0 + 28.0 * sqrt21) / 45.0,
     -(273.0 + 53.0 * sqrt21) / 72.0,
     (301.0 + 53.0 * sqrt21) / 72.0,
     (28.0 - 28.0 * sqrt21) / 45.0,
     (49.0 - 7.0 * sqrt21) / 18.0}};

/* Prepare the workspace for the integration. */
void *
verner_alloc(size_t number_of_equations)
{
    verner_state_t *state = (verner_state_t *)malloc(sizeof(verner_state_t));
    size_t i;

    state->y_temp = (double *)malloc(number_of_equations * sizeof(double));

    for (i = 0; i < NUMBER_OF_STAGES; i++)
        state->k[i] = (double *)malloc(number_of_equations * sizeof(double));

    return state;
}

void
verner_compute_stages(void *vstate,
                      size_t number_of_equations,
                      double t,
                      double h_step,
                      double y[number_of_equations],
                      void (*RHS)(double, double *, double *))
{
    verner_state_t *state = (verner_state_t *)vstate;

    size_t i;

    double *const y_temp = state->y_temp;

    /* Obtain the memory locations of the stages k1 through k11.
     Note that k1 is stored in state->k[0] due to zero-based indexing. */
    double *const k1 = state->k[0];
    double *const k2 = state->k[1];
    double *const k3 = state->k[2];
    double *const k4 = state->k[3];
    double *const k5 = state->k[4];
    double *const k6 = state->k[5];
    double *const k7 = state->k[6];
    double *const k8 = state->k[7];
    double *const k9 = state->k[8];
    double *const k10 = state->k[9];
    double *const k11 = state->k[10];

    /* k1 stage */
    RHS(t, y, k1);

    /* k2 stage */
    for (i = 0; i < number_of_equations; i++)
        y_temp[i] = y[i] + h_step * (RK_A[1][0] * k1[i]);
    RHS(t + h_step * RK_c[1], y_temp, k2);

    /* k3 stage */
    for (i = 0; i < number_of_equations; i++)
        y_temp[i] = y[i] + h_step * (RK_A[2][0] * k1[i] + RK_A[2][1] * k2[i]);
    RHS(t + h_step * RK_c[2], y_temp, k3);

    /* k4 stage */
    for (i = 0; i < number_of_equations; i++)
        y_temp[i] =
            y[i] + h_step * (RK_A[3][0] * k1[i] + RK_A[3][1] * k2[i]
            + RK_A[3][2] * k3[i]);
    RHS(t + h_step * RK_c[3], y_temp, k4);

    /* k5 stage */
    for (i = 0; i < number_of_equations; i++)
        y_temp[i] =
            y[i] + h_step * (RK_A[4][0] * k1[i] + RK_A[4][1] * k2[i]
            + RK_A[4][2] * k3[i] + RK_A[4][3] * k4[i]);
    RHS(t + h_step * RK_c[4], y_temp, k5);

    /* k6 stage */
    for (i = 0; i < number_of_equations; i++)
        y_temp[i] =
            y[i] + h_step * (RK_A[5][0] * k1[i] + RK_A[5][1] * k2[i]
            + RK_A[5][2] * k3[i] + RK_A[5][3] * k4[i] + RK_A[5][4] * k5[i]);
    RHS(t + h_step * RK_c[5], y_temp, k6);

    /* k7 stage */
    for (i = 0; i < number_of_equations; i++)
        y_temp[i] =
            y[i] + h_step * (RK_A[6][0] * k1[i] + RK_A[6][1] * k2[i]
            + RK_A[6][2] * k3[i] + RK_A[6][3] * k4[i] + RK_A[6][4] * k5[i]
            + RK_A[6][5] * k6[i]);
    RHS(t + h_step * RK_c[6], y_temp, k7);

    /* k8 stage */
    for (i = 0; i < number_of_equations; i++)
        y_temp[i] =
            y[i] + h_step * (RK_A[7][0] * k1[i] + RK_A[7][1] * k2[i]
            + RK_A[7][2] * k3[i] + RK_A[7][3] * k4[i] + RK_A[7][4] * k5[i]
            + RK_A[7][5] * k6[i] + RK_A[7][6] * k7[i]);
    RHS(t + h_step * RK_c[7], y_temp, k8);

    /* k9 stage */
    for (i = 0; i < number_of_equations; i++)
        y_temp[i] =
            y[i] + h_step * (RK_A[8][0] * k1[i] + RK_A[8][1] * k2[i]
            + RK_A[8][2] * k3[i] + RK_A[8][3] * k4[i] + RK_A[8][4] * k5[i]
            + RK_A[8][5] * k6[i] + RK_A[8][6] * k7[i] + RK_A[8][7] * k8[i]);
    RHS(t + h_step * RK_c[8], y_temp, k9);

    /* k10 stage */
    for (i = 0; i < number_of_equations; i++)
        y_temp[i] =
            y[i] + h_step * (RK_A[9][0] * k1[i] + RK_A[9][1] * k2[i]
            + RK_A[9][2] * k3[i] + RK_A[9][3] * k4[i] + RK_A[9][4] * k5[i]
            + RK_A[9][5] * k6[i] + RK_A[9][6] * k7[i] + RK_A[9][7] * k8[i]
            + RK_A[9][8] * k9[i]);
    RHS(t + h_step * RK_c[9], y_temp, k10);

    /* k11 stage */
    for (i = 0; i < number_of_equations; i++)
        y_temp[i] =
            y[i] + h_step * (RK_A[10][0] * k1[i] + RK_A[10][1] * k2[i]
            + RK_A[10][2] * k3[i] + RK_A[10][3] * k4[i] + RK_A[10][4] * k5[i]
            + RK_A[10][5] * k6[i] + RK_A[10][6] * k7[i] + RK_A[10][7] * k8[i]
            + RK_A[10][8] * k9[i] + RK_A[10][9] * k10[i]);
    RHS(t + h_step * RK_c[10], y_temp, k11);

    return;
}

void
verner_dense_output(void *vstate,
                    size_t number_of_equations,
                    double t,
                    double theta,
                    double h_step,
                    double y[number_of_equations])
{
    verner_state_t *state = (verner_state_t *)vstate;

    size_t i;

    /* Fetch the stages k1 and k8 through k11. */
    double *const k1 = state->k[0];
    double *const k8 = state->k[7];
    double *const k9 = state->k[8];
    double *const k10 = state->k[9];
    double *const k11 = state->k[10];

    static const double interp_coeffs[5][5] = {
        {1.0, 
        -5.0, 
         10.0, 
        -35.0 / 4.0,
        14.0 / 5.0},
        {0.0, 
        (49.0 + 7.0 * sqrt21) / 12.0, 
        (-245.0 - 21.0 * sqrt21) / 18.0, 
        (196.0 + 7.0 * sqrt21) / 12.0, 
        -98.0 / 15.0},
        {0.0,
        -8.0 / 3.0,
        128.0 / 9.0,
        -56.0 / 3.0,
        112.0 / 15.0},
        {0.0, 
        (49.0 - 7.0 * sqrt21) / 12.0, 
        (-245.0 + 21.0 * sqrt21) / 18.0, 
        (196.0 - 7.0 * sqrt21) / 12.0, 
        -98.0 / 15.0},
        {0.0, 
        -1.0 / 2.0, 
        3.0, 
        -21.0 / 4.0, 
        14.0 / 5.0}
    };

    /* Calculate the basis functions required for the interpolation. */
    double b_function[5];
    b_function[0] = interp_coeffs[0][0] + theta * (interp_coeffs[0][1]
    + theta * (interp_coeffs[0][2] + theta * (interp_coeffs[0][3]
    + theta * (interp_coeffs[0][4]))));
    b_function[1] = theta * (interp_coeffs[1][1] + theta * (interp_coeffs[1][2]
    + theta * (interp_coeffs[1][3] + theta * (interp_coeffs[1][4]))));
    b_function[2] = theta * (interp_coeffs[2][1] + theta * (interp_coeffs[2][2]
    + theta * (interp_coeffs[2][3] + theta * (interp_coeffs[2][4]))));
    b_function[3] = theta * (interp_coeffs[3][1] + theta * (interp_coeffs[3][2]
    + theta * (interp_coeffs[3][3] + theta * (interp_coeffs[3][4]))));
    b_function[4] = theta * (interp_coeffs[4][1] + theta * (interp_coeffs[4][2]
    + theta * (interp_coeffs[4][3] + theta * (interp_coeffs[4][4]))));

    /* Combine the k-stages in the final sum to approximate the solution. */
    double y_output[number_of_equations];
    for (i = 0; i < number_of_equations; i++)
        y_output[i] = y[i] + (theta * h_step) * (b_function[0] * k1[i]
        + b_function[1] * k8[i] + b_function[2] * k9[i] + b_function[3] * k10[i]
        + b_function[4] * k11[i]);

    /* Print the results to console. */
    printf("t = %8.6lf", t + (theta * h_step));
    for (i = 0; i < number_of_equations; i++)
        printf(", y[%li] = %20.16lf", i, y_output[i]);

    printf("\n");

    return;
}

void
verner_apply(void *vstate,
             size_t number_of_equations, double h_step, double y[number_of_equations])
{
    verner_state_t *state = (verner_state_t *)vstate;

    double *const k1 = state->k[0];
    double *const k8 = state->k[7];
    double *const k9 = state->k[8];
    double *const k10 = state->k[9];
    double *const k11 = state->k[10];

    /* Combine the k-stages in the final sum to forward the solution. */
    for (size_t i = 0; i < number_of_equations; i++)
        y[i] = y[i] + h_step * (RK_b[0] * k1[i] + RK_b[7] * k8[i] + RK_b[8] * k9[i]
        + RK_b[9] * k10[i] + RK_b[10] * k11[i]);

    return;
}

void
verner_free(void *vstate)
{
    verner_state_t *state = (verner_state_t *)vstate;

    for (size_t i = 0; i < NUMBER_OF_STAGES; i++)
        free(state->k[i]);

    free(state->y_temp);
    free(state);
}
