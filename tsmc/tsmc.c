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
    Hmm hmm;
    Hmm_init(&hmm, n, ts);
    Hmm_set_lambdas(&hmm, n, lambdas);
    Hmm_make_omega_intervals(&hmm);
    Hmm_get_pis(&hmm);
    Hmm_get_expectations(&hmm);
    Hmm_get_qts(&hmm);
    Hmm_get_pts(&hmm);
    Hmm_get_emissions(&hmm);
    //Hmm_print_emissions(&hmm);

    Data dat;
    Data_init(&dat, polarized);

    FILE * fin = chfopen("testseqs", "r");
    Data_read_data(&dat, fin);
    fclose(fin);

    Em em;
    Em_init(&em, &dat, &hmm);
    Em_get_forward(&em);
    Em_get_backward(&em);
    Em_get_expectations(&em);

    //Data_print_seqs(&dat);
    
    Em_free(&em);
    Data_free(&dat);
    Hmm_free(&hmm);
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
