#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "hmm.h"
#include "definitions.h"
#include "data.h"
#include "em.h"
#include "options.h"

Em em;

int main(int argc, char ** argv)
{

    int i;
    double F;
    const int n = 10;
    double * ts = (double *)chmalloc(sizeof(double) * (n+1));
    ts[0] = 0.0;
    double * lambdas = (double *)chmalloc(sizeof(double) * (n+1));
    lambdas[0] = 1.0;
    for(i = 1; i < n+1; i++)
    {
        F = (double)(i) * 1.0/((double)n+1.0);
        ts[i] = -log(1-F);
        lambdas[i] = (i <= 3) ? 1.0 : 2.0;
    }

    double initTheta = 1e-3;
    double initRho = 1e-3;
    double initTd = 0.2;

    Hmm hmm;
    Hmm_init(&hmm, n, ts);
    Hmm_make_hmm(&hmm, lambdas, n, initTheta, initRho, initTd);
    Hmm_print_qts(&hmm);
    chfree(ts);
    chfree(lambdas);

    //Hmm_print_pis(&hmm);
    //Hmm_print_expectations(&hmm);
    //Hmm_print_pts(&hmm);
    //Hmm_free(&hmm);
    return 0;
}
