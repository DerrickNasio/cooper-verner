/*
 * This routine performs dense output for the Cooper-Verner algorithm.
 * A 4th-order sum of basis functions is used as the interpolant.
 */
#include "Verner.h"

void DenseOutput(int n_equat, double t, double y[n_equat], double h_step, double t_output[], struct workspace *myWorkspace)
{
    /* Fetch the stages k1 and k8 through k11. */
    double *const k1 = myWorkspace->k[0];
    double *const k8 = myWorkspace->k[7];
    double *const k9 = myWorkspace->k[8];
    double *const k10 = myWorkspace->k[9];
    double *const k11 = myWorkspace->k[10];

    int outputSize = *(&t_output + 1) - t_output;

    for (int i = 0; i < outputSize; i++)
    {
        /* Calculate the auxiliary variable theta such that 0 <= theta <= 1 */
        double theta = (t_output[i] - t) / h_step;

        double basisFunctions[5];
        if (theta >= 0.0 && theta <= 1.0)
        {
            basisFunctions[0] = 1.0 + theta * (-5.0 + theta * (10.0 + theta * (-35.0 / 4.0 + theta * (14.0 / 5.0))));
            basisFunctions[1] = theta * ((49.0 + 7.0 * sqrt21) / 12.0 + theta * ((-245.0 - 21.0 * sqrt21) / 18.0 + theta * ((196.0 + 7.0 * sqrt21) / 12.0 + theta * (-98.0 / 15.0))));
            basisFunctions[2] = theta * (-8.0 / 3.0 + theta * (128.0 / 9.0 + theta * (-56.0 / 3.0 + theta * (112.0 / 15.0))));
            basisFunctions[3] = theta * ((49.0 - 7.0 * sqrt21) / 12.0 + theta * ((-245.0 + 21.0 * sqrt21) / 18.0 + theta * ((196.0 - 7.0 * sqrt21) / 12.0 + theta * (-98.0 / 15.0))));
            basisFunctions[4] = theta * (-1.0 / 2.0 + theta * (3.0 + theta * (-21.0 / 4.0 + theta * (14.0 / 5.0))));

            double y_output[n_equat];
            /* Combine the k-stages in the final sum to approximate the solution. */
            for (int j = 0; j < n_equat; j++)
                y_output[j] = y[j] + (theta * h_step) * (basisFunctions[0] * k1[j] + basisFunctions[1] * k8[j] + basisFunctions[2] * k9[j] + basisFunctions[3] * k10[j] + basisFunctions[4] * k11[j]);

            /* Print the results to console. */
            printf("t = %8.lf", t + (theta * h_step));
            for (int j = 0; j < n_equat; j++)
                printf(", y[%i] = %20.16lf", j, y_output[j]);

            printf("\n");
        }
        else
        {
            continue;
        }
    }
    return;
}
