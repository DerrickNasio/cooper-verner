/*
 * This routine performs dense output for the Cooper-Verner algorithm.
 * A 4th-order sum of basis functions is used as the interpolant.
 */
#include "Verner.h"

void dense_output(
    size_t number_of_equations,
    double t,
    double theta,
    double h_step,
    double y[number_of_equations],
    struct workspace *workspace_struct)
{
    /* Fetch the stages k1 and k8 through k11. */
    double *const k1 = workspace_struct->k[0];
    double *const k8 = workspace_struct->k[7];
    double *const k9 = workspace_struct->k[8];
    double *const k10 = workspace_struct->k[9];
    double *const k11 = workspace_struct->k[10];

    /* Calculate the five basis functions required for the interpolation. */
    double b_function[5];
    b_function[0] = 1.0 + theta * (-5.0 + theta * (10.0 + theta * (-35.0 / 4.0 + theta * (14.0 / 5.0))));
    b_function[1] = theta * ((49.0 + 7.0 * sqrt21) / 12.0 + theta * ((-245.0 - 21.0 * sqrt21) / 18.0 + theta * ((196.0 + 7.0 * sqrt21) / 12.0 + theta * (-98.0 / 15.0))));
    b_function[2] = theta * (-8.0 / 3.0 + theta * (128.0 / 9.0 + theta * (-56.0 / 3.0 + theta * (112.0 / 15.0))));
    b_function[3] = theta * ((49.0 - 7.0 * sqrt21) / 12.0 + theta * ((-245.0 + 21.0 * sqrt21) / 18.0 + theta * ((196.0 - 7.0 * sqrt21) / 12.0 + theta * (-98.0 / 15.0))));
    b_function[4] = theta * (-1.0 / 2.0 + theta * (3.0 + theta * (-21.0 / 4.0 + theta * (14.0 / 5.0))));

    /* Combine the k-stages in the final sum to approximate the solution. */
    double y_output[number_of_equations];
    for (int j = 0; j < number_of_equations; j++)
        y_output[j] = y[j] + (theta * h_step) * (b_function[0] * k1[j] + b_function[1] * k8[j] + b_function[2] * k9[j] + b_function[3] * k10[j] + b_function[4] * k11[j]);

    /* Print the results to console. */
    printf("t = %8.lf", t + (theta * h_step));
    for (int j = 0; j < number_of_equations; j++)
        printf(", y[%i] = %20.16lf", j, y_output[j]);

    printf("\n");

    return;
}
