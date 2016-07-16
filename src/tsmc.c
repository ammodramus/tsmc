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
    get_ts_psmc(ts, opt.maxT, n);

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
    timestamp("reading data");
    Data_read_data(&dat, fin);
    timestamp("data read");
    if(strcmp(opt.filename, "-") != 0)
    {
        fclose(fin);
    }

    double initRho;
    timestamp("getting initial rho...");
    initRho = Em_get_initial_rho(&dat);
    timestamp("finished initial rho, initializing em...");

    Em_init(&em, &dat, ts, initRho, opt.numFreeLambdas, opt.n,
            opt.numEmIterations, opt.numOptimizations, opt.lambdaCounts,
            opt.asexEnabled, opt.flagDt);
    timestamp("em initialized");

    for(i = 0; i < opt.numEmIterations; i++)
    {
        timestamp("starting iteration");
        Em_iterate(&em);
        timestamp("ending iteration");
        Em_print_iteration(&em);
    }

    Em_free(&em);
    Em_free_globals();
    Data_free(&dat);
    chfree(ts);

    return 0;
}
