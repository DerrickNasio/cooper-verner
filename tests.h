/*
 * tests.h
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// tests.c
void RHSLinear(double t, double y[], double dydt[]);
void RHSExp(double t, double y[], double dydt[]);
void RHSSineCosine(double t, double y[], double dydt[]);
void RHSCosine(double t, double y[], double dydt[]);
void RHSStepFunction(double t, double y[], double dydt[]);
void RHSClassicStiff(double t, double y[], double dydt[]);
void RHSvanderPol(double t, double y[], double dydt[]);
void RHSStiffTrig(double t, double y[], double dydt[]);
