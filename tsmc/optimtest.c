#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "definitions.h"
#include "asa047.h"

#define N 10

// minimized at f(1,1,1,1,1,1,1,...) = 0
double rosenbrock(double * par)
{
    double sum = 0.0;
    int i;
    for(i = 0; i < N-1; i++)
    {
        sum += 100 * (par[i+1]-par[i]*par[i])*(par[i+1]-par[i]*par[i]) +
            (par[i]-1.0)*(par[i]-1.0);
    }
    return sum;
}

int main(int argc, char ** argv)
{
    const int numParams = N;
    int i;

    double * start = (double *)chmalloc(sizeof(double) * numParams);
    double * end = (double *)chmalloc(sizeof(double) * numParams);
    double * step = (double *)chmalloc(sizeof(double) * numParams);
    for(i = 0; i < N; i++)
    {
        start[i] = 0.0;
        step[i] = 0.1;
    }
    double fmin;

    int konvge = 1, maxNumEval = 1000000;
    int iterationCount = 0, numRestarts, errorNum;
    double reqmin = 1e-16;

    nelmin(rosenbrock, numParams, start, end,
        &fmin, reqmin, step, konvge, maxNumEval, &iterationCount,
        &numRestarts, &errorNum);

    for(i = 0; i < N; i++)
    {
        printf("%i: %g\n", i, end[i]);
    }
    printf("fmin: %g\n", fmin);
    return 0;
}
