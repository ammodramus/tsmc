#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "hmm.h"
#include "definitions.h"
#include "data.h"
#include "em.h"

Em em;

int main(int argc, char ** argv)
{
    int i, n = 10;
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

    Data dat;
    Data_init(&dat, polarized);
    FILE * fin;
    if(argc == 2)
    {
        fin = chfopen(argv[1], "r");
    }
    else
    {
        fin = chfopen("testseqs", "r");
    }
    Data_read_data(&dat, fin);
    fclose(fin);

    Em_init(&em, &dat, lambdas, ts, n, initTheta, initRho, initTd, maxIterations);

    for(i = 0; i < 2; i++)
    {
        Em_iterate(&em);
    }

    // Em object should have two Hmm objects as members, one current, one that
    // is maximized.  Have a flag that switches which is which, iteratively.
    // each iteration:
    //  - get forward
    //  - get backward
    //  - get expectations
    //  - maximize likelihood
    //  - flip flag

    //Data_print_seqs(&dat);
    
    Em_free(&em);
    Data_free(&dat);
    free(ts);
    free(lambdas);

    //Hmm_print_demography(&hmm);
    //Hmm_print_pis(&hmm);
    //Hmm_print_expectations(&hmm);
    //Hmm_print_qts(&hmm);
    //Hmm_print_pts(&hmm);
    //Hmm_free(&hmm);
    return 0;
}
