#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "definitions.h"
#include "genotypelikelihoods.h"
#include "hmm.h"
#include "definitions.h"

int main()
{
    int polarized, bin_width, n, error, i;
    double maxT, initTd;
    FILE *fin;
    polarized = 0;
    GenotypeLikeData dat;
    GenotypeLikeData_init(&dat, polarized);

    polarized = 0;
    fin = fopen("test_likelihoods.txt", "r");
    GenotypeLikeData_read_data(fin, &dat, 0);

    n = 15;
    maxT = 2.0;
    initTd = 0.3;

    double * ts = (double *)chmalloc(sizeof(double) * (n+1));
    get_ts_psmc(ts, maxT, n);


    double * lambdas = (double *)chmalloc(sizeof(double) * (n+1));
    for(i = 0; i < n+1; i++)
    {
        lambdas[i] = 1.0;
    }

    Hmm hmm;
    Hmm_init(&hmm, n);
    Hmm_make_hmm(&hmm, lambdas, ts, n, 0.2, 0.2, initTd, &error);

    bin_width = 10;
    GenLikeEmissions emissions;
    GenLikeEmissions_calculate(&emissions, &dat, bin_width, &hmm);
    GenLikeEmissions_print_first(&emissions);

    GenLikeEmissions_free(&emissions);
    GenotypeLikeData_free(&dat);
    Hmm_free(&hmm);
    chfree(ts);
    chfree(lambdas);
    fclose(fin);

    return 0;
}
