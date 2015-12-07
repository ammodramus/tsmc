#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define DEBUGREPORTI(x) fprintf(stderr, #x " = %i\n", x)
#define DEBUGREPORTF(x) fprintf(stderr, #x " = %f\n", x)

typedef struct 
{
    int n; // number of changepoints
    int numStates;
    double * ts; // 
    double * lambdas;
    double * deltas;
    double * intervalOmegas;
    double * pis;
    double * Eijs3s;
    double * Eijs2s;
    double ** qts;
} Hmm;

void * chmalloc(size_t size)
{
    void * ptr = malloc(size);
    if(!ptr)
    {
        fprintf(stderr, "Out of memory.\n");
        exit(1);
    }
    return ptr;
}

inline int get_index(int i, int j, int n)
{
    assert(i>=0 && j>=0 && n>=0 && i <= n && j <= n && i <= j);
    //int idx = i*n-i*(i-1)/2+j+1;
    const int idx = i*n-i*(i-1)/2+j;
    assert(idx < (n+1)*(n+2)/2);
    return idx;
}

void Hmm_init(Hmm * hmm, int n, double * ts)
{
    int i;
    const int numStates = (n+1)*(n+2)/2;

    assert(n > 0);

    hmm->n = n;
    hmm->lambdas = (double *)chmalloc(sizeof(double) * (n+1));
    hmm->deltas = (double *)chmalloc(sizeof(double) * n);
    hmm->intervalOmegas = (double *)chmalloc(sizeof(double) * n);
    hmm->ts = (double *)chmalloc(sizeof(double) * (n+1));
    hmm->Eijs2s = (double *)chmalloc(sizeof(double)*numStates);
    hmm->Eijs3s = (double *)chmalloc(sizeof(double)*numStates);
    hmm->pis = (double *)chmalloc(sizeof(double)*numStates);
    hmm->qts = (double **)chmalloc(sizeof(double *)*numStates);
    hmm->numStates = numStates;
    for(i = 0; i < numStates; i++)
    {
        hmm->qts[i] = (double *)chmalloc(sizeof(double)*numStates);
        hmm->Eijs3s[i] = -1.0;
        hmm->Eijs2s[i] = -1.0;
        hmm->pis[i] = -1.0;
    }
    
    for(i = 0; i < n+1; i++)
    {
        hmm->lambdas[i] = 1.0;
        hmm->ts[i] = ts[i];
    }
    for(i = 0; i < n; i++)
    {
        hmm->deltas[i] = hmm->ts[i+1]-hmm->ts[i];
    }
    return;
}

void Hmm_free(Hmm * hmm)
{
    int i;
    free(hmm->lambdas);
    free(hmm->deltas);
    free(hmm->ts);
    free(hmm->Eijs3s);
    free(hmm->Eijs2s);
    free(hmm->intervalOmegas);
    free(hmm->pis);
    for(i = 0; i < hmm->numStates; i++)
    {
        free(hmm->qts[i]);
    }
    free(hmm->qts);
    return;
}


void Hmm_set_lambdas(Hmm * hmm, const int n, const double * lambdas)
{
    int i;
    assert(n == hmm->n);
    for(i = 0; i < n+1; i++)
    {
        hmm->lambdas[i] = lambdas[i];
    }
    return;
}

/* 
 * which to make first
 * \Omega coefficients
 * \pi's
 * \Eij's
 * q(t|s)'s
 */

void Hmm_make_omega_intervals(Hmm * hmm)
{
    int i;
    const int n = hmm->n;
    for(i = 0; i < n; i++)
    {
        hmm->intervalOmegas[i] = hmm->deltas[i]/hmm->lambdas[i];
    }
    return;
}

double get_omega_interval_interval(Hmm * hmm, int a, int b)
{
    int i;
    double omega = 0.0;
    assert(a <= b);
    if(a == b)
    {
        return 0.0;
    }
    for(i = a; i < b; i++)
    {
        omega += hmm->intervalOmegas[i];
    }

    assert(omega >= 0);

    return omega;
}

// Calculate the equilibrium joint probabilities of s_3 and s_2 being in the 
// state (i,j)
void Hmm_get_pis(Hmm * hmm)
{
    const int n = hmm->n;
    int i, j, idx;
    for(i = 0; i < n+1; i++)
    {
        for(j = i; j < n+1; j++)
        {
            idx = get_index(i, j, n);
            if(i == j)
            {
                if(i == n)
                {
                    hmm->pis[idx] = exp(-3.0*get_omega_interval_interval(
                                hmm, 0, n));
                }
                else
                {
                    hmm->pis[idx] = 1.0/2.0*exp(-3.0*get_omega_interval_interval(
                            hmm, 0, i))*
                            (2.0 - 3.0*exp(-hmm->intervalOmegas[i]) 
                             + exp(-3.0*hmm->intervalOmegas[i]));
                }
            }
            // i < j
            else
            {
                assert(i<j);
                if(j < n)
                {
                    hmm->pis[idx] = 3.0/2.0*
                        exp(-3.0*get_omega_interval_interval(hmm,0,i))
                        * exp(-get_omega_interval_interval(hmm,i+1,j))
                        * (exp(-hmm->intervalOmegas[i])
                                - exp(-3.0*hmm->intervalOmegas[i]))
                        * (1.0 - exp(-hmm->intervalOmegas[j]));
                }
                else
                {
                    assert(j == n);
                    hmm->pis[idx] = 3.0/2.0*
                        exp(-3.0*get_omega_interval_interval(hmm,0,i)) *
                        exp(-get_omega_interval_interval(hmm,i+1,j)) *
                        (exp(-hmm->intervalOmegas[i]) - 
                                exp(-3.0*hmm->intervalOmegas[i]));
                }
            }
        }
    }
    return;
}

void Hmm_get_expectations(Hmm * hmm)
{
    const int n = hmm->n;
    int i, j, idx;
    for(i = 0; i < n+1; i++)
    {
        for(j = i; j < n+1; j++)
        {
            idx = get_index(i, j, n);
            if(i == j)
            {
                if(i == n)
                {
                    hmm->Eijs3s[idx] = hmm->ts[n] + hmm->lambdas[n]/3.0;
                    hmm->Eijs2s[idx] = hmm->Eijs3s[idx] + hmm->lambdas[n];
                }
                else
                {
                    assert(i < n && i == j);
                    hmm->Eijs3s[idx] = 1.0/(12.0*hmm->pis[idx])
                        * exp(-3*get_omega_interval_interval(hmm,0,i))
                        * (4*(3*hmm->ts[i] + hmm->lambdas[i]) +
                        exp(-3*hmm->intervalOmegas[i]) *
                         (6*hmm->ts[i+1]+5*hmm->lambdas[i]) -
                        9*exp(-hmm->intervalOmegas[i]) *
                        (2*hmm->ts[i] + hmm->lambdas[i]));
                        
                    hmm->Eijs2s[idx] = 1.0/(6.0*hmm->pis[idx])
                        * exp(-3*get_omega_interval_interval(hmm, 0, i))
                        * (exp(-3*hmm->intervalOmegas[i])
                                * (3*hmm->ts[i+1]+hmm->lambdas[i])
                           + 6*hmm->ts[i] + 8*hmm->lambdas[i]
                           - 9*exp(-hmm->intervalOmegas[i]) *
                           (hmm->ts[i+1]+hmm->lambdas[i]));
                }
            }
            else
            {
                assert(i < j);
                if(j == n)
                {
                    hmm->Eijs3s[idx] = 3.0/(4.0*hmm->pis[idx]) *
                      exp(-3*get_omega_interval_interval(hmm,0,i))*
                      exp(-get_omega_interval_interval(hmm,i+1,n))*
                      ((hmm->lambdas[i]+2*hmm->ts[i])*
                       exp(-hmm->intervalOmegas[i]) -
                       (hmm->lambdas[i]+2*hmm->ts[i+1])*
                       exp(-3*hmm->intervalOmegas[i]));
                    
                    hmm->Eijs2s[idx] = 3.0/(2.0*hmm->pis[idx]) *
                        exp(-3*get_omega_interval_interval(hmm,0,i))*
                        exp(-get_omega_interval_interval(hmm,i+1,n))*
                        (exp(-hmm->intervalOmegas[i]) -
                         exp(-3*hmm->intervalOmegas[i]))*
                        (hmm->lambdas[n] + hmm->ts[n]);
                }
                else
                {
                    assert(j < n && i < j);
                    hmm->Eijs3s[idx] = 3.0/(4.0*hmm->pis[idx])*
                        exp(-3*get_omega_interval_interval(hmm,0,i))*
                        (1.0-exp(-hmm->intervalOmegas[j]))*
                        exp(-get_omega_interval_interval(hmm,i+1,j))*
                        ((hmm->lambdas[i]+2*hmm->ts[i])*
                         exp(-hmm->intervalOmegas[i]) -
                         (hmm->lambdas[i]+2*hmm->ts[i+1])*
                         exp(-3*hmm->intervalOmegas[i]));
                    hmm->Eijs2s[idx] = 3.0/(2.0*hmm->pis[idx])*
                        exp(-3*get_omega_interval_interval(hmm,0,i))*
                        exp(-get_omega_interval_interval(hmm,i+1,j))*
                        (exp(-hmm->intervalOmegas[i])-
                         exp(-3*hmm->intervalOmegas[i]))*
                        (hmm->lambdas[j]+hmm->ts[j]-
                         (hmm->lambdas[j]+hmm->ts[j+1])*
                         exp(-hmm->intervalOmegas[j]));
                }
            }
        }
    }
    return;
}

// calculate \Omega(\E_{Es3_i, Es3_j}[s_3], T_b) 
double get_omega_Es3_interval(Hmm * hmm, const int b, const int Es3_i, 
        const int Es3_j)
{
    int i;
    double omega = 0.0;

    assert(Es3_i < b);
    assert(b >= 0 && Es3_i >= 0 && Es3_j >= 0);
    // all less than n because we are taking interval above
    assert(b < hmm->n && Es3_i < hmm->n && Es3_j < hmm->n);

    const Es3idx = get_index(Es3_i, Es3_j, hmm->n);
    const double Es3 = hmm->Eijs3s[Es3idx];

    omega += (hmm->ts[Es3_i+1] - Es3) * hmm->intervalOmegas[Es3_i];
    omega += get_omega_interval_interval(hmm, Es3_i+1, b);

    assert(omega >= 0);

    return omega;
}

// calculate \Omega(\E_{Es2_i, Es2_j}[s_2], T_b) 
double get_omega_Es2_interval(Hmm * hmm, const int b, const int Es2_i, 
        const int Es2_j)
{
    int i;
    double omega = 0.0;

    assert(Es2_j < b);
    assert(Es2_i <= Es2_j);
    assert(Es2_i >= 0 && b >= 0 && Es2_i >= 0 && Es2_j >= 0);
    // all less than n because we are taking interval above
    assert(Es2_i < hmm->n && b < hmm->n && Es2_j < hmm->n);

    const Es2idx = get_index(Es2_i, Es2_j, hmm->n);
    const double Es2 = hmm->Eijs2s[Es2idx];

    omega += (hmm->ts[Es2_j+1] - Es2) * hmm->intervalOmegas[Es2_j];
    omega += get_omega_interval_interval(hmm, Es2_j+1, b);

    assert(omega >= 0);

    return omega;
}

double get_omega_interval_Es3(Hmm * hmm, const int a, const int Es3_i, 
        const int Es3_j)
{
    int i;
    double omega = 0.0;
    assert(a < Es3_i); // maybe we want <=
    assert(Es3_i <= Es3_j);
    assert(Es3_i >= 0 && a >= 0 && Es3_i >= 0 && Es3_j >= 0);
    assert(Es3_i <= hmm->n && a < hmm->n && Es3_j <= hmm->n);

    const Es3idx = get_index(Es3_i, Es3_j, hmm->n);
    const double Es3 = hmm->Eijs3s[Es3idx];

    omega += get_omega_interval_interval(hmm, a, Es3_i);
    omega += (Es3 - hmm->ts[Es3_i])/hmm->lambdas[Es3_i];

    assert(omega >= 0);

    return omega;
}

double get_omega_interval_Es2(Hmm * hmm, const int a, const int Es2_i,
        const int Es2_j)
{
    int i;
    double omega = 0.0;

    assert(a < Es2_j); // maybe we want <=
    assert(Es2_i <= Es2_j);
    assert(a >= 0 && Es2_i >= 0 && Es2_j >= 0);
    assert(Es2_i <= hmm->n && a < hmm->n && Es2_j <= hmm->n);

    const Es2idx = get_index(Es2_i, Es2_j, hmm->n);
    const double Es2 = hmm->Eijs2s[Es2idx];

    omega += get_omega_interval_interval(hmm, a, Es2_j);
    omega += (Es2 - hmm->ts[Es2_j])/hmm->lambdas[Es2_j];

    assert(omega >= 0);

    return omega;
}

double get_omega_Es3_Es2(Hmm * hmm, const int Es3_i, const int Es3_j)
{
    int i;
    double omega = 0.0;

    assert(Es3_i <= Es3_j);
    assert(Es3_i >= 0 && Es3_j >= 0);
    assert(Es3_i <= hmm->n && Es3_j <= hmm->n);

    const Es3idx = get_index(Es3_i, Es3_j, hmm->n);
    const double Es3 = hmm->Eijs3s[Es3idx];

    if(Es3_i == Es3_j)
    {
        return (hmm->Eijs2s[Es3idx] - hmm->Eijs3s[Es3idx])/hmm->lambdas[Es3_i];
    }

    omega += get_omega_Es3_interval(hmm, Es3_j, Es3_i, Es3_j);
    omega += (Es3 - hmm->ts[Es3_j])/hmm->lambdas[Es3_j];

    assert(omega >= 0);
    return omega;
}

void Hmm_print_demography(Hmm * hmm)
{
    int i;
    for(i = 0; i < hmm->n+1; i++)
        printf("%i\t%f\t%f\n", i, hmm->ts[i], hmm->lambdas[i]);
    return;
}

void Hmm_print_pis(Hmm * hmm)
{
    int i, j, idx;
    const int n = hmm->n;
    for(i = 0; i < n+1; i++)
    {
        for(j = i; j < n+1; j++)
        {
            idx = get_index(i, j, n);
            printf("%i\t%i\t%f\t%f\t%f\t%f\t%f\n",
                i, j, hmm->ts[i],
                (i < n) ? hmm->ts[i+1] : INFINITY, hmm->ts[j], 
                (j < n) ? hmm->ts[j+1] : INFINITY, hmm->pis[idx]);
        }
    }
    return;
}

void Hmm_print_expectations(Hmm * hmm)
{
    int i, j, idx;
    const int n = hmm->n;
    for(i = 0; i < n+1; i++)
    {
        for(j = i; j < n+1; j++)
        {
            idx = get_index(i, j, n);
            printf("%i\t%i\t%f\t%f\t%f\t%f\t%f\t%f\n",
                i, j, hmm->ts[i],
                (i < n) ? hmm->ts[i+1] : INFINITY, hmm->ts[j], 
                (j < n) ? hmm->ts[j+1] : INFINITY,
                hmm->Eijs3s[idx], hmm->Eijs2s[idx]);
        }
    }
    return;
}

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
    Hmm_print_demography(&hmm);
    //Hmm_print_pis(&hmm);
    //Hmm_print_expectations(&hmm);
    //Hmm_free(&hmm);
    return 0;
}
