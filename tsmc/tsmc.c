#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "hmm.h"
#include "definitions.h"
#include "data.h"
#include "em.h"
#include "options.h"
#include "random.h"

Em em;

int main(int argc, char ** argv)
{

    Options opt;
    Options_parse_options(&opt, argc, argv);

    randseed();

    const int n = opt.n;

    int i;
    double * ts = (double *)chmalloc(sizeof(double) * (n+1));
    get_ts(ts, n);

    Data dat;
    Data_init(&dat, polarized);
    FILE * fin;
    // strcmp returns 0 with match
    if(strcmp(opt.filename, "-") != 0)
    {
        fin = chfopen(opt.filename, "r");
    }
    else
    {
        fin = stdin;
    }
    Data_read_data(&dat, fin);
    if(strcmp(opt.filename, "-") != 0)
    {
        fclose(fin);
    }

    double initRho;
    initRho = Em_get_initial_rho(&dat);
    const double initTd = 0.2;

    Em_init(&em, &dat, ts, initRho, initTd, opt.numFreeLambdas, opt.n,
            opt.numEmIterations, opt.lambdaCounts);

    for(i = 0; i < opt.numEmIterations; i++)
    {
        Em_iterate(&em);
        Em_print_iteration(&em);
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

    //Hmm_print_demography(&hmm);
    //Hmm_print_pis(&hmm);
    //Hmm_print_expectations(&hmm);
    //Hmm_print_qts(&hmm);
    //Hmm_print_pts(&hmm);
    //Hmm_free(&hmm);
    return 0;
}
