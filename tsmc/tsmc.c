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
    if(opt.psmcIntervals)
    {
        get_ts_psmc(ts, opt.maxT, n);
    }
    else
    {
        get_ts_msmc(ts, n);
    }

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
    timestamp("getting initial rho...");
    initRho = Em_get_initial_rho(&dat);
    timestamp("finished initial rho...");
    const double asexInitTd = 0.2;
    const double initTd = opt.asexEnabled ? asexInitTd : 0.0;

    Em_init(&em, &dat, ts, initRho, initTd, opt.numFreeLambdas, opt.n,
            opt.numEmIterations, opt.lambdaCounts, opt.asexEnabled);

    for(i = 0; i < opt.numEmIterations; i++)
    {
        Em_iterate(&em);
        Em_print_iteration(&em);
    }

    Em_free(&em);
    Data_free(&dat);
    free(ts);

    return 0;
}
