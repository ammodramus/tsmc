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
    FILE * fin = chfopen("diptriptestseqs", "r");
    Data_read_data(&dat, fin);

    double initRho = 0.04;

    Em_init(&em, &dat, ts, initRho, opt.numFreeLambdas, opt.n,
            opt.numEmIterations, opt.numOptimizations, opt.lambdaCounts,
            1, 0);
    Em_iterate(&em);
    Em_print_iteration(&em);
    Em_free(&em);

    Em_init(&em, &dat, ts, initRho, opt.numFreeLambdas, opt.n,
            opt.numEmIterations, opt.numOptimizations, opt.lambdaCounts,
            1, 1);
    for(i = 0; i < 10; i++)
    {
        Em_iterate(&em);
        Em_print_iteration(&em);
    }

    double D3, lamd;
    double D3min = -1.0 * em.hmm[em.hmmFlag].Td;
    double D3max = 0.0;
    double lamdMin = 1e-10;
    double lamdMax = 1e3;

    double D3step = (D3max - D3min)/(100.0);

    double expect;

    double * par = (double *)chmalloc(sizeof(double) * 6);
    par[0] = sqrt(em.hmm[em.hmmFlag].rho);
    par[1] = sqrt(em.hmm[em.hmmFlag].theta);
    par[2] = sqrt(em.hmm[em.hmmFlag].maxT);

    //const double rho = par[0]*par[0];
    //const double theta = par[1]*par[1];
    //const double maxT = par[2]*par[2];
    //const double dipTime = par[3]*par[3]; // *amount* of time in diploid asex state
    //const double tripTime = par[4]*par[4]; // *amount* of time in triploid asex state
    //const double lamd = par[5]*par[5];
    //
    em.hmmFlag = !em.hmmFlag;

    const int numNonLamParams = 6;
    const int numParams = em.numFreeLambdas + numNonLamParams;
    for(i = numNonLamParams; i < numParams; i++)
    {
        par[i] = sqrt(em.freeLambdas[i-numNonLamParams]);
    }

    for(D3 = D3min; D3 <= D3max; D3 += D3step)
    {
        for(lamd = lamdMin; lamd <= lamdMax; lamd *= 10.0)
        {
            par[3] = sqrt(-1.0*D3);
            par[4] = sqrt(em.hmm[em.hmmFlag].Td + D3);
            par[5] = sqrt(lamd);
            expect = objective_function_dt(par);
            printf("%.15f\t%.15f\t%.15f\n", D3, lamd, expect);
        }
    }

    chfree(par);
    Em_free(&em);
    Data_free(&dat);
    chfree(ts);

    return 0;
}
