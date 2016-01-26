#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "hmm.h"
#include "definitions.h"
#include "data.h"
#include "em.h"

int main(int argc, char ** argv)
{
    const int n = 10;
    int i;
    double F;
    double * ts = (double *)chmalloc(sizeof(double) * (n+1));
    double * lambdas = (double *)chmalloc(sizeof(double) * (n+1));
    
    ts[0] = 0.0;
    lambdas[0] = 1.0;
    for(i = 1; i < n+1; i++)
    {
        F = (double)(i) * 1.0/((double)n+1.0);
        ts[i] = -log(1-F);
        if(i < 4)
        {
            lambdas[i] = 1.0;
        }
        else
        {
            lambdas[i] = 2.0;
        }
    }

    double initRho = 1e-3;
    double initTheta = 1e-3;
    double initTd = 0.2;
    const int maxIterations = 20;

    Hmm hmm;
    Hmm_init(&hmm, n);
    int error;
    Hmm_make_hmm(&hmm, lambdas, ts, n, initTheta, initRho, initTd, &error);
    Hmm_print_pis(&hmm);

    return 0;
}
