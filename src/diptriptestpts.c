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
    double theta = 0.4;
    double rho = 0.2;
    double Td = 0.2;
    double D3 = 0.0;
    double lamd = 1.0;

    Hmm hmm;
    Hmm_init(&hmm, n);
    Hmm hmmdt;
    Hmm_init_dt(&hmmdt, n);

    int error;
    Hmm_make_hmm(&hmm, lambdas, ts, n, theta, rho, Td, &error);

    int errorDt;
    Hmm_make_hmm_dt(&hmmdt, lambdas, ts, n, theta, rho, Td, D3, lamd, &error);

    Hmm_print_pts(&hmm);
    printf("---\n");
    Hmm_print_pts(&hmmdt);

    return 0;
}
