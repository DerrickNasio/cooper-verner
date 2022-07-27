/*
 * tests.c
 * Some test functions to be used with the integrator.
 */
#include "tests.h"

/* y'(t) = 2.
 * Solution: y(t) = y(0) + 2 * t.
 */
void RHSLinear(double t, double y[], double dydt[])
{
	dydt[0] = 2.0;
}

/* y'(t) = y(t).
 * Initial value: y(0) = 1.0.
 * Solution: y(t) = exp(t)
 */
void RHSExp(double t, double y[], double dydt[])
{
	dydt[0] = y[0];
}

/* y0'(t) = - y1(t), y1'(t) = y0(t).
 * Initial values: y0(0) = 1.0, y1(0) = 0.0.
 * Solution: y0(t) = cos(t), y1(t) = sin(t).
 */
void RHSSineCosine(double t, double y[], double dydt[])
{
	dydt[0] = -y[1];
	dydt[1] = y[0];
}

/* y'(t) = cos(t).
 * Initial value: y(0) = 0.0.
 * Solution: y(t) = sin(t)
 */
void RHSCosine(double t, double y[], double dydt[])
{
	dydt[0] = cos(t);
}

/* RHS for step function:
 * Derivative change at t = 1.0
 * from y'(t) = 0.0 to y'(t) = 1.0.
 */
void RHSStepFunction(double t, double y[], double dydt[])
{
	if (t >= 1.0)
		dydt[0] = 1;
	else
		dydt[0] = 0;
}

/* RHS for classic stiff example:
 * y0'(t) =  998 * y0(t) + 1998 * y1(t),
 * y1'(t) = -999 * y0(t) - 1999 * y1(t).
 * Initial values: y0(0) = 1.0, y1(0) = 0.0.
 * Solution:
 * y0(t) = 2 * exp(-t) - exp(-1000 * t),
 * y1(t) = - exp(-t) + exp(-1000 * t).
 */
void RHSClassicStiff(double t, double y[], double dydt[])
{
	dydt[0] = 998.0 * y[0] + 1998.0 * y[1];
	dydt[1] = -999.0 * y[0] - 1999.0 * y[1];
}

/* RHS for van der Pol oscillator:
 * y0'(t) = y1(t),
 * y1'(t) = - y0(t) + mu * y1(t) * (1 - y0(t)^2).
 * Initial values: y0(0) = 1.0, y1(0) = 0.0.
 */
void RHSvanderPol(double t, double y[], double dydt[])
{
	const double mu = 10.0;

	dydt[0] = y[1];
	dydt[1] = -y[0] + mu * y[1] * (1.0 - y[0] * y[0]);
}

/* RHS for stiff trigonometric example:
 * y'(t) = -50 * (y(t) - cos(t)).
 * Initial value: y0(0) = 0.0.
 */
void RHSStiffTrig(double t, double y[], double dydt[])
{
	dydt[0] = -50 * (y[0] - cos(t));
}
