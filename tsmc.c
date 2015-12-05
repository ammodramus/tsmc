#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

typedef struct 
{
    int n; // number of changepoints
    double * ts; // 
    double * lambdas;
    double * deltas;
    double * intervalOmegas;
    double * Eijs3s;
    double * Eijs2s;
} Hmm;

inline int index(int i, int j, int n)
{
    return i*n-i*(i-1)/2+j+1;
}

void Hmm_init(Hmm * hmm, int n, double * ts)
{
    int i;
    const int numStates = (n+1)*(n+2)/2;

    assert(n > 0);

    hmm->n = n;
    // TODO replace these mallocs with a custom, checking malloc
    hmm->lambdas = (double *)malloc(sizeof(double) * (n+1));
    hmm->deltas = (double *)malloc(sizeof(double) * n);
    hmm->intervalOmegas = (double *)malloc(sizeof(double) * n);
    hmm->ts = (double *)malloc(sizeof(double) * (n+1));
    hmm->Eijs2s = (double *)malloc(sizeof(double)*numStates);
    hmm->Eijs3s = (double *)malloc(sizeof(double)*numStates);
    
    
    for(i = 0; i < n+1; i++)
    {
        hmm->lambdas[i] = 1.0;
        hmm->ts[i] = ts[i];
    }
    for(i = 0; i < n; i++)
        hmm->deltas[i] = hmm->ts[i+1]-hmm->ts[i];
    return;
}

void Hmm_set_lambdas(Hmm * hmm, const int n, const double * lambdas)
{
    int i;
    assert(n == hmm->n);
    for(i = 0; i < n+1; i++)
        hmm->lambdas[i] = lambdas[i];
    return;
}

/* 
 * which to make first
 * \Omega coefficients
 * \pi's
 * \Eij's
 * q(t|s)'s
 */
void Hmm_make_hmm(Hmm * hmm)
{
    return;
}


void Hmm_print_demography(Hmm * hmm)
{
    int i;
    for(i = 0; i < hmm->n+1; i++)
        printf("%i\t%f\t%f\n", i, hmm->ts[i], hmm->lambdas[i]);
    return;
}

void Hmm_free(Hmm * hmm)
{
    free(hmm->lambdas);
    free(hmm->deltas);
    free(hmm->ts);
    free(hmm->Eijs3s);
    free(hmm->Eijs2s);
    return;
}

// a list of quantities that are computed again and again
// e^{-k\Omega(0,T_i)}
// e^{-k\Omega(T_i,T_j)}, i<j
// e^{-\Delta_i/\lambda_i}
// \E_{i,j}[s_3], \E_{i,j}[s_2]
// make a table of the cumulative coalescent rates between i and j

int main(int argc, char ** argv)
{
    int i, n = 10;
    double F;
    double * ts = (double *)malloc(sizeof(double) * (n+1));
    
    ts[0] = 0.0;
    for(i = 1; i < n+1; i++)
    {
        F = (double)(i) * 1.0/((double)n+1.0);
        ts[i] = -log(1-F);
    }
    Hmm hmm;
    Hmm_init(&hmm, n, ts);
    Hmm_print_demography(&hmm);
    Hmm_free(&hmm);
    return 0;
}
