#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "definitions.h"
#include "hmm.h"
#include "qtscases.h"
#include "qtssuppdt.h"

inline int get_index(const int i, const int j, const int n)
{
    assert(0 <= i && i <= j && j <= n);
    const int idx = i*n-i*(i-1)/2+j;
    assert(idx < (n+1)*(n+2)/2);
    return idx;
}

inline int get_index_dt(const int i, const int j, const int n, const int W)
{
    assert(0 <= i && i <= j && j <= n && (W == 0 || W == 1));
    const int idx = i*n-i*(i-1)/2+j + W*(n+1)*(n+2)/2;
    assert(0 <= idx && idx < (n+1)*(n+2));
    return idx;
}

void Hmm_init(Hmm * hmm, const int n)
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
    hmm->pts = (double **)chmalloc(sizeof(double *)*numStates);
    hmm->emissions = (sixd *)chmalloc(sizeof(sixd)*numStates); 
    hmm->numStates = numStates;
    for(i = 0; i < numStates; i++)
    {
        hmm->qts[i] = (double *)chmalloc(sizeof(double)*numStates);
        hmm->pts[i] = (double *)chmalloc(sizeof(double)*numStates);
        hmm->Eijs3s[i] = -1.0;
        hmm->Eijs2s[i] = -1.0;
        hmm->pis[i] = -1.0;
    }
    hmm->theta = -1.0;
    hmm->rho = -1.0;
    hmm->Td = -1.0;
    hmm->maxT = -1.0;

    hmm->flagDt = 0;    // set flag
    hmm->numStatesDt = -1;
    hmm->lambdaDt = -1.0;
    hmm->D3 = -1.0;

    hmm->qtsDt = NULL;

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
        free(hmm->pts[i]);
    }
    free(hmm->qts);
    free(hmm->pts);
    free(hmm->emissions);
    return;
}

void Hmm_init_dt(Hmm * hmm, const int n)
{
    int i;
    const int numStates = (n+1)*(n+2)/2;
    const int numStatesDt = 2*numStates;

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
    hmm->qtsDt = (double **)chmalloc(sizeof(double *)*numStatesDt);
    hmm->pts = (double **)chmalloc(sizeof(double *)*numStatesDt);
    hmm->emissions = (sixd *)chmalloc(sizeof(sixd)*numStatesDt); 
    hmm->numStates = numStates;
    hmm->numStatesDt = numStatesDt;
    for(i = 0; i < numStatesDt; i++)
    {
        hmm->qtsDt[i] = (double *)chmalloc(sizeof(double)*numStatesDt);
        hmm->pts[i] = (double *)chmalloc(sizeof(double)*numStatesDt);
    }
    for(i = 0; i < numStates; i++)
    {
        hmm->qts[i] = (double *)chmalloc(sizeof(double)*numStates);
        hmm->Eijs3s[i] = -1.0;
        hmm->Eijs2s[i] = -1.0;
        hmm->pis[i] = -1.0;
    }
    hmm->theta = -1.0;
    hmm->rho = -1.0;
    hmm->Td = -1.0;
    hmm->maxT = -1.0;
    hmm->lambdaDt = -1.0;
    hmm->D3 = -1.0;

    hmm->flagDt = 1;    // set flag

    return;
}

void Hmm_free_dt(Hmm * hmm)
{
    int i;
    free(hmm->lambdas);
    free(hmm->deltas);
    free(hmm->ts);
    free(hmm->Eijs3s);
    free(hmm->Eijs2s);
    free(hmm->intervalOmegas);
    free(hmm->pis);
    for(i = 0; i < hmm->numStatesDt; i++)
    {
        free(hmm->qtsDt[i]);
        free(hmm->pts[i]);
    }
    for(i = 0; i < hmm->numStates; i++)
    {
        free(hmm->qts[i]);
    }
    free(hmm->qts);
    free(hmm->qtsDt);
    free(hmm->pts);
    free(hmm->emissions);
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

void Hmm_set_ts_and_deltas(Hmm * hmm, const double * ts)
{
    const int n = hmm->n;
    int i;
    for(i = 0; i < n+1; i++)
    {
        hmm->ts[i] = ts[i];
    }
    for(i = 0; i < n; i++)
    {
        hmm->deltas[i] = hmm->ts[i+1]-hmm->ts[i];
    }
    hmm->maxT = ts[n];
    return;
}

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

double get_omega_interval_interval(Hmm * hmm, const int a, const int b)
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
                    assert(0 < hmm->pis[idx] && hmm->pis[idx] < 1);
                }
                else
                {
                    hmm->pis[idx] = 1.0/2.0*exp(-3.0*get_omega_interval_interval(
                            hmm, 0, i))*
                            (2.0 - 3.0*exp(-hmm->intervalOmegas[i]) 
                             + exp(-3.0*hmm->intervalOmegas[i]));
                    assert(0 < hmm->pis[idx] && hmm->pis[idx] < 1);
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
                    assert(0 < hmm->pis[idx] && hmm->pis[idx] < 1);
                }
                else
                {
                    assert(j == n);
                    hmm->pis[idx] = 3.0/2.0*
                        exp(-3.0*get_omega_interval_interval(hmm,0,i)) *
                        exp(-get_omega_interval_interval(hmm,i+1,j)) *
                        (exp(-hmm->intervalOmegas[i]) - 
                                exp(-3.0*hmm->intervalOmegas[i]));
                    assert(0 < hmm->pis[idx] && hmm->pis[idx] < 1);
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
                    assert(0 < hmm->Eijs3s[idx] && hmm->Eijs3s[idx] < INFINITY);
                    assert(0 < hmm->Eijs2s[idx] && hmm->Eijs2s[idx] < INFINITY);
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
                    assert(0 < hmm->Eijs3s[idx] && hmm->Eijs3s[idx] < INFINITY);
                        
                    hmm->Eijs2s[idx] = 1.0/(6.0*hmm->pis[idx])
                        * exp(-3*get_omega_interval_interval(hmm, 0, i))
                        * (exp(-3*hmm->intervalOmegas[i])
                                * (3*hmm->ts[i+1]+hmm->lambdas[i])
                           + 6*hmm->ts[i] + 8*hmm->lambdas[i]
                           - 9*exp(-hmm->intervalOmegas[i]) *
                           (hmm->ts[i+1]+hmm->lambdas[i]));
                    assert(0 < hmm->Eijs2s[idx] && hmm->Eijs3s[idx] < INFINITY);
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
    assert(b >= 0 && Es3_i >= 0 && Es3_j >= Es3_i);
    assert(b <= hmm->n && Es3_i < hmm->n && Es3_j <= hmm->n);

    const Es3idx = get_index(Es3_i, Es3_j, hmm->n);
    const double Es3 = hmm->Eijs3s[Es3idx];

    omega += (hmm->ts[Es3_i+1] - Es3) / hmm->lambdas[Es3_i];
    omega += get_omega_interval_interval(hmm, Es3_i+1, b);

    assert(omega >= 0);

    return omega;
}

// calculate \Omega(\E_{Es2_i, Es2_j}[s_2], T_b) 
double get_omega_Es2_interval(Hmm * hmm, const int b, const int Es2_i, 
        const int Es2_j)
{
    assert(0 <= Es2_i && Es2_i <= Es2_j && Es2_j < b && b <= hmm->n);

    const Es2idx = get_index(Es2_i, Es2_j, hmm->n);
    const double Es2 = hmm->Eijs2s[Es2idx];

    double omega = 0.0;
    omega += (hmm->ts[Es2_j+1] - Es2) / hmm->lambdas[Es2_j];
    omega += get_omega_interval_interval(hmm, Es2_j+1, b);

    assert(omega >= 0);

    return omega;
}

double get_omega_interval_Es3(Hmm * hmm, const int a, const int Es3_i, 
        const int Es3_j)
{
    int i;
    double omega = 0.0;
    assert(0 <= a && a <= Es3_i && Es3_i <= Es3_j && Es3_j <= hmm->n);

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

    assert(0 <= a && a <= Es2_j && Es2_i <= Es2_j && Es2_j <= hmm->n);

    const Es2idx = get_index(Es2_i, Es2_j, hmm->n);
    const double Es2 = hmm->Eijs2s[Es2idx];

    omega += get_omega_interval_interval(hmm, a, Es2_j);
    omega += (Es2 - hmm->ts[Es2_j])/hmm->lambdas[Es2_j];

    assert(omega >= 0);

    return omega;
}

double get_omega_Es3_Es2(Hmm * hmm, const int Es3_i, const int Es3_j)
{
    double omega = 0.0;

    assert(0 <= Es3_i && Es3_i <= Es3_j && Es3_j <= hmm->n);

    const Es3idx = get_index(Es3_i, Es3_j, hmm->n);
    const double Es3 = hmm->Eijs3s[Es3idx];
    const double Es2 = hmm->Eijs2s[Es3idx];

    if(Es3_i == Es3_j)
    {
        omega = (hmm->Eijs2s[Es3idx] - hmm->Eijs3s[Es3idx]) / hmm->lambdas[Es3_i];
        assert(omega >= 0);
        return omega;
    }

    omega += get_omega_Es3_interval(hmm, Es3_j, Es3_i, Es3_j);
    omega += (Es2 - hmm->ts[Es3_j])/hmm->lambdas[Es3_j];

    assert(omega >= 0);
    return omega;
}

void Hmm_get_qts(Hmm * hmm)
{
    int i, j, k, l, rowIdx, colIdx;
    const int n = hmm->n;
    const int numStates = hmm->numStates;

    double ** const qts = hmm->qts;

    for(i=0; i<n+1; i++)
    {
        for(j=i; j<n+1; j++)
        {
            for(k=0; k<n+1; k++)
            {
                for(l=k; l<n+1; l++)
                {
                    rowIdx = get_index(i,j,n);
                    colIdx = get_index(k,l,n);
                    qts[rowIdx][colIdx] = 0.0; // as a default
                    if(i == k && k < j && j < l)
                    {
                        qts[rowIdx][colIdx] = qts_case_A(hmm, i, j, k, l);
                    }
                    else if(i == k && k < l && l < j)
                    {
                        qts[rowIdx][colIdx] = qts_case_B(hmm, i, j, k, l);
                    }
                    else if(k < i && i == l && l < j)
                    {
                        qts[rowIdx][colIdx] = qts_case_C(hmm, i, j, k, l);
                    }
                    else if(k < i && i < j && j == l)
                    {
                        qts[rowIdx][colIdx] = qts_case_D(hmm, i, j, k, l);
                    }
                    else if(i < k && k < j && j == l)
                    {
                        qts[rowIdx][colIdx] = qts_case_E(hmm, i, j, k, l);
                    }
                    else if(i < j && j == k && k < l)
                    {
                        qts[rowIdx][colIdx] = qts_case_F(hmm, i, j, k, l);
                    }
                    else if(i == k && k == l && l < j)
                    {
                        qts[rowIdx][colIdx] = qts_case_G(hmm, i, j, k, l);
                        qts[rowIdx][colIdx] += qts_case_G2(hmm, i, j, k, l);
                    }
                    else if(i < k && k == l && l == j)
                    {
                        qts[rowIdx][colIdx] = qts_case_H(hmm, i, j, k, l);
                        qts[rowIdx][colIdx] += qts_case_H2(hmm, i, j, k, l);
                    }
                    else if(i == j && j == k && k < l)
                    {
                        qts[rowIdx][colIdx] = qts_case_I(hmm, i, j, k, l);
                        qts[rowIdx][colIdx] += qts_case_I2(hmm, i, j, k, l);
                    }
                    else if(k < i && i == j && j == l)
                    {
                        qts[rowIdx][colIdx] = qts_case_J(hmm, i, j, k, l);
                        qts[rowIdx][colIdx] += qts_case_J2(hmm, i, j, k, l);
                    }
                }
            }
        }
    }
    return;
}

void Hmm_get_pts(Hmm * hmm)
{
    const int numStates = hmm->numStates;
    const int n = hmm->n;
    const double rho = hmm->rho;
    const double theta = hmm->theta;
    double ** const qts = hmm->qts;
    double ** const pts = hmm->pts;
    double * const Es3s = hmm->Eijs3s;
    double * const Es2s = hmm->Eijs2s;

    assert(!hmm->flagDt);

    int i, j, k, l, ijidx, klidx;
    double Es3, Es2, treeSize, probRecomb, sum;

    for(i = 0; i <= n; i++)
    {
        for(j = i; j <= n; j++)
        {
            ijidx = get_index(i, j, n);
            Es3 = Es3s[ijidx];
            Es2 = Es2s[ijidx];
            treeSize = 3.0*Es3 + 2.0*(Es2-Es3);
            probRecomb = 1.0-exp(-treeSize * rho/2.0);

            sum = 0.0;
            
            for(k = 0; k <= n; k++)
            {
                for(l = k; l <=n; l++)
                {
                    klidx = get_index(k, l, n);
                    if(klidx == ijidx)
                    {
                        continue;
                    }
                    pts[ijidx][klidx] = qts[ijidx][klidx] * probRecomb;
                    sum += pts[ijidx][klidx];
                }
            }
            pts[ijidx][ijidx] = 1.0 - sum;
        }
    }

    return;
}

void Hmm_get_qts_dt(Hmm * hmm)
{
    int W, Wp, i, j, k, l, rowIdx, colIdx, expIdx;
    const int n = hmm->n;
    const int numStates = hmm->numStates;

    double ** const qts = hmm->qtsDt;

    double Es3, Es2;

    double ratio;

    assert(hmm->flagDt);

    //Hmm_get_qts(hmm);
    
    for(W = 0; W <= 1; W++)
    {
        for(Wp = 0; Wp <= 1; Wp++)
        {
            for(i=0; i<n+1; i++)
            {
                for(j=i; j<n+1; j++)
                {
                    rowIdx = get_index_dt(i,j,n,W);
                    expIdx = get_index(i,j,n);
                    Es3 = hmm->Eijs3s[expIdx];
                    Es2 = hmm->Eijs2s[expIdx];
                    ratio = (2*Es2 + Es3)/(2*Es2+Es3 - hmm->D3);
                    for(k=0; k<n+1; k++)
                    {
                        for(l=k; l<n+1; l++)
                        {
                            colIdx = get_index_dt(k,l,n,Wp);
                            qts[rowIdx][colIdx] = 0.0; // as a default
                            if(i == k && k < j && j < l) // case A
                            {
                                if(W == Wp)
                                {
                                    qts[rowIdx][colIdx] = ratio * qts_case_A(hmm, i, j, k, l);
                                    if(W == 0)
                                    {
                                        assert(Wp == 0);
                                        qts[rowIdx][colIdx] += qts_case_A_supp(hmm, i, j, k, l);
                                    }
                                }
                            }
                            else if(i == k && k < l && l < j) // case B
                            {
                                if(W == Wp)
                                {
                                    qts[rowIdx][colIdx] = ratio * qts_case_B(hmm, i, j, k, l);
                                    if(W == 0)
                                    {
                                        assert(Wp == 0);
                                        qts[rowIdx][colIdx] += qts_case_B_supp(hmm, i, j, k, l);
                                    }
                                }
                            }
                            else if(k < i && i == l && l < j) // case C
                            {
                                if(W == 0 && Wp == 1)
                                {
                                    qts[rowIdx][colIdx] = ratio * qts_case_C(hmm, i, j, k, l) +
                                        qts_case_C_supp(hmm, i, j, k, l);
                                }
                                if(W == 1)
                                {
                                    qts[rowIdx][colIdx] = 1.0/2.0 * ratio * qts_case_C(hmm, i, j, k, l);
                                }
                            }
                            else if(k < i && i < j && j == l) // case D
                            {
                                if(W == Wp)
                                {
                                    qts[rowIdx][colIdx] = ratio * qts_case_D(hmm, i, j, k, l);
                                    if(W == 1)
                                    {
                                        assert(Wp == 1);
                                        qts[rowIdx][colIdx] += qts_case_D_supp(hmm, i, j, k, l);
                                    }
                                }
                            }
                            else if(i < k && k < j && j == l) // case E
                            {
                                if(W == Wp)
                                {
                                    qts[rowIdx][colIdx] = ratio * qts_case_E(hmm, i, j, k, l);
                                    if(W == 1)
                                    {
                                        assert(Wp == 1);
                                        qts[rowIdx][colIdx] += qts_case_E_supp(hmm, i, j, k, l);
                                    }
                                }
                            }
                            else if(i < j && j == k && k < l) // case F
                            {
                                if(W == 1)
                                {
                                    qts[rowIdx][colIdx] = 1.0/2.0 * ratio * qts_case_F(hmm, i, j, k, l);
                                    if(Wp == 0)
                                    {
                                        qts[rowIdx][colIdx] += qts_case_F_supp(hmm, i, j, k, l);
                                    }
                                }
                                if(W == 0 && Wp == 1)
                                {
                                    qts[rowIdx][colIdx] = ratio * qts_case_F(hmm, i, j, k, l);
                                }
                            }
                            else if(i == k && k == l && l < j) // case G and G2
                            {
                                qts[rowIdx][colIdx] = 0.0;
                                if(W == Wp) // case G
                                {
                                    qts[rowIdx][colIdx] += ratio * qts_case_G(hmm, i, j, k, l);
                                    if(W == 0)
                                    {
                                        assert(Wp == 0);
                                        qts[rowIdx][colIdx] += qts_case_G_supp(hmm, i, j, k, l);
                                    }
                                }
                                if(W == 0 && Wp == 1) // case G2
                                {
                                    qts[rowIdx][colIdx] += ratio * qts_case_G2(hmm,i,j,k,l) +
                                        qts_case_G2_supp(hmm,i,j,k,l);
                                }
                                if(W == 1) // also case G2
                                {
                                    qts[rowIdx][colIdx] += 0.5 * ratio * qts_case_G2(hmm, i,j,k,l);
                                }
                            }
                            else if(i < k && k == l && l == j) // case H and H2
                            {
                                qts[rowIdx][colIdx] = 0.0;
                                if(W == Wp) // case H
                                {
                                    qts[rowIdx][colIdx] += ratio * qts_case_H(hmm,i,j,k,l);
                                    if(W == 1)
                                    {
                                        assert(Wp == 1);
                                        qts[rowIdx][colIdx] += qts_case_H_supp(hmm,i,j,k,l);
                                    }

                                }
                                if(W == 0 && Wp == 1) // case H2
                                {
                                    qts[rowIdx][colIdx] += ratio * qts_case_H2(hmm,i,j,k,l);
                                }
                                if(W == 1) // also case H2
                                {
                                    qts[rowIdx][colIdx] += 0.5 * ratio * qts_case_H2(hmm,i,j,k,l);
                                    if(Wp == 0)
                                    {
                                        qts[rowIdx][colIdx] += qts_case_H2_supp(hmm,i,j,k,l);
                                    }
                                }
                            }
                            else if(i == j && j == k && k < l) // case I and I2
                            {
                                qts[rowIdx][colIdx] = 0.0;
                                if(W == Wp) // case I
                                {
                                    qts[rowIdx][colIdx] += ratio * qts_case_I(hmm,i,j,k,l);
                                    if(W == 0)
                                    {
                                        assert(Wp == 0);
                                        qts[rowIdx][colIdx] += qts_case_I_supp(hmm,i,j,k,l);
                                    }
                                }
                                if(W == 0 && Wp == 1) // case I2
                                {
                                    qts[rowIdx][colIdx] += ratio * qts_case_I2(hmm,i,j,k,l);
                                }
                                if(W == 1) // also case I2
                                {
                                    qts[rowIdx][colIdx] += 0.5*ratio*qts_case_I2(hmm,i,j,k,l);
                                    if(Wp == 0)
                                    {
                                        qts[rowIdx][colIdx] += qts_case_I2_supp(hmm,i,j,k,l);
                                    }
                                }
                            }

                            else if(k < i && i == j && j == l) // cases J and J2
                            {
                                qts[rowIdx][colIdx] = 0.0;
                                if(W == Wp) // case J
                                {
                                    qts[rowIdx][colIdx] += ratio * qts_case_J(hmm,i,j,k,l);
                                    if(W == 1)
                                    {
                                        assert(Wp == 1);
                                        qts[rowIdx][colIdx] += qts_case_J_supp(hmm,i,j,k,l);
                                    }
                                }
                                if(W == 0 && Wp == 1)
                                {
                                    qts[rowIdx][colIdx] += ratio * qts_case_J2(hmm,i,j,k,l) +
                                        qts_case_J2_supp(hmm,i,j,k,l);
                                }
                                if(W == 1)
                                {
                                    qts[rowIdx][colIdx] += 0.5 * ratio * qts_case_J2(hmm,i,j,k,l);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return;
}

void Hmm_get_pts_dt(Hmm * hmm)
{
    const int n = hmm->n;
    const double rho = hmm->rho;
    const double theta = hmm->theta;
    const double D3 = hmm->D3;
    double ** const qts = hmm->qtsDt;
    double ** const pts = hmm->pts;
    double * const Es3s = hmm->Eijs3s;
    double * const Es2s = hmm->Eijs2s;


    assert(hmm->flagDt);

    int W, Wp, i, j, k, l, ijidx, klidx, expIdx;
    double Es3, Es2, treeSize, probRecomb, sum;

    for(W = 0; W <= 1; W++)
    {
        for(i = 0; i <= n; i++)
        {
            for(j = i; j <= n; j++)
            {
                ijidx = get_index_dt(i, j, n, W);
                expIdx = get_index(i, j, n);
                Es3 = Es3s[expIdx];
                Es2 = Es2s[expIdx];
                treeSize = 3.0*Es3 + 2.0*(Es2-Es3) - D3;
                probRecomb = 1.0-exp(-treeSize * rho/2.0);

                sum = 0.0;
                
                for(Wp = 0; Wp <= 1; Wp++)
                {
                    for(k = 0; k <= n; k++)
                    {
                        for(l = k; l <=n; l++)
                        {
                            klidx = get_index_dt(k, l, n, Wp);
                            if(klidx == ijidx)
                            {
                                continue;
                            }
                            pts[ijidx][klidx] = qts[ijidx][klidx] * probRecomb;
                            sum += pts[ijidx][klidx];
                        }
                    }
                    assert(0.0 <= sum && sum <= 1.0);
                    pts[ijidx][ijidx] = 1.0 - sum;
                }
            }
        }
    }
    return;
}

void Hmm_get_emissions(Hmm * hmm)
{
    const double numStates = hmm->numStates;
    const double n = hmm->n;
    const double theta = hmm->theta;
    const double Td = hmm->Td;
    double * const Es3s = hmm->Eijs3s;
    double * const Es2s = hmm->Eijs2s;
    sixd * const emissions = hmm->emissions;

    assert(!hmm->flagDt);

    int i, j, ijidx;
    double Es3, Es2, treeSize;

    for(i = 0; i <= n; i++)
    {
        for(j = i; j <= n; j++)
        {
            ijidx = get_index(i,j,n);
            Es3 = Es3s[ijidx];
            Es2 = Es2s[ijidx];

            treeSize = 2.0*Es2 + Es3;

            /* homozygous */
            emissions[ijidx][0] = exp(-(treeSize + 3.0*Td)*theta/2.0);
            /* 1-mut */
            emissions[ijidx][1] = exp(-(Es2-Es3)*theta/2.0) *
                (1 - exp(-(2*Es3 + Es2 + 3*Td)*theta/2.0));
            /* 2-mut */
            emissions[ijidx][2] = (1.0 - exp(-(Es2-Es3)*theta/2.0)) * 
                exp(-(2*Es3 + Es2 + 3*Td)*theta/2.0);
            /* 1-mut and 2-mut */
            emissions[ijidx][3] = (1.0 - exp(-(Es2-Es3)*theta/2.0)) * 
                (1 - exp(-(2*Es3 + Es2 + 3*Td)*theta/2.0));
            /* heterozygous, no polarization */
            emissions[ijidx][4] = 1.0 - emissions[ijidx][0];
            /* missing */
            emissions[ijidx][5] = 1.0;
        }
    }

    return;
}

void Hmm_get_emissions_dt(Hmm * hmm)
{
    const double numStates = hmm->numStates;
    const double n = hmm->n;
    const double theta = hmm->theta;
    const double Td = hmm->Td;
    double * const Es3s = hmm->Eijs3s;
    double * const Es2s = hmm->Eijs2s;
    sixd * const emissions = hmm->emissions;

    assert(hmm->flagDt);

    int W, i, j, ijidx, expIdx;
    double Es3, Es2, treeSize;

    for(W = 0; W <= 1; W++)
    {
        for(i = 0; i <= n; i++)
        {
            for(j = i; j <= n; j++)
            {
                ijidx = get_index_dt(i,j,n,W);
                expIdx = get_index(i,j,n);
                Es3 = Es3s[expIdx];
                Es2 = Es2s[expIdx];

                treeSize = 2.0*Es2 + Es3;

                emissions[ijidx][0] = exp(-(treeSize + 3.0*Td)*theta/2.0);
                emissions[ijidx][1] = exp(-(Es2-Es3)*theta/2.0) *
                    (1 - exp(-(2*Es3 + Es2 + 3*Td)*theta/2.0));
                emissions[ijidx][2] = (1.0 - exp(-(Es2-Es3)*theta/2.0)) * 
                    exp(-(2*Es3 + Es2 + 3*Td)*theta/2.0);
                emissions[ijidx][3] = (1.0 - exp(-(Es2-Es3)*theta/2.0)) * 
                    (1 - exp(-(2*Es3 + Es2 + 3*Td)*theta/2.0));
                emissions[ijidx][4] = 1.0 - emissions[ijidx][0];
                emissions[ijidx][5] = 1.0;
            }
        }
    }

    return;
}

inline void Hmm_set_theta(Hmm * hmm, double theta)
{
    hmm->theta = theta;
    return;
}

inline void Hmm_set_rho(Hmm * hmm, double rho)
{
    hmm->rho = rho;
    return;
}

inline void Hmm_set_Td(Hmm * hmm, double Td)
{
    hmm->Td = Td;
    return;
}

inline void Hmm_set_D3(Hmm * hmm, double D3)
{
    assert(hmm->flagDt);
    hmm->D3 = D3;
}

inline void Hmm_set_lambdaD(Hmm * hmm, double lamd)
{
    assert(hmm->flagDt);
    hmm->lambdaDt = lamd;
    return;
}

void Hmm_make_hmm(Hmm * hmm, double * lambdas, double * ts,
        int numChangepoints, double theta, double rho, double Td, int * error)
{

    *error = 0;

    // numChangepoints is n in the notation of the paper
    assert(hmm->n == numChangepoints);
    Hmm_set_ts_and_deltas(hmm, ts);
    Hmm_set_lambdas(hmm, numChangepoints, lambdas);
    Hmm_set_theta(hmm, theta);
    Hmm_set_rho(hmm, rho);
    Hmm_set_Td(hmm, Td);

    Hmm_make_omega_intervals(hmm);
    int i;
    for(i = 0; i < hmm->n; i++)
    {
        if(hmm->intervalOmegas[i] > 10)
        {
            *error = 1;
            return;
        }
    }

    Hmm_get_pis(hmm);
    Hmm_get_expectations(hmm);
    Hmm_get_qts(hmm);
    Hmm_get_pts(hmm);
    Hmm_get_emissions(hmm);

    return;
}

void Hmm_make_hmm_dt(Hmm * hmm, double * lambdas, double * ts,
        int numChangepoints, double theta, double rho, double Td, double D3,
        int * error)
{

    *error = 0;

    // numChangepoints is n in the notation of the paper
    assert(hmm->n == numChangepoints);
    Hmm_set_ts_and_deltas(hmm, ts);
    Hmm_set_lambdas(hmm, numChangepoints, lambdas);
    Hmm_set_theta(hmm, theta);
    Hmm_set_rho(hmm, rho);
    Hmm_set_Td(hmm, Td);
    Hmm_set_D3(hmm, D3);
    Hmm_set_lambdaD(hmm, 1.0);

    Hmm_make_omega_intervals(hmm);
    int i;
    for(i = 0; i < hmm->n; i++)
    {
        if(hmm->intervalOmegas[i] > 10)
        {
            *error = 1;
            return;
        }
    }

    Hmm_get_pis(hmm);
    Hmm_get_expectations(hmm);
    Hmm_get_qts_dt(hmm);
    Hmm_get_pts_dt(hmm);
    Hmm_get_emissions_dt(hmm);

    return;
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

void Hmm_print_qts(Hmm * hmm)
{
    assert(hmm);
    int row, col;
    const int n = hmm->n;
    const int numStates = hmm->numStates;
    for(row = 0; row < numStates; row++)
    {
        for(col = 0; col < numStates-1; col++)
        {
            printf("%f,", hmm->qts[row][col]);
        }
        printf("%f\n", hmm->qts[row][numStates-1]);
    }
    return;
}

void Hmm_print_pts(Hmm * hmm)
{
    assert(hmm);
    int row, col;
    int numStates;
    if(!hmm->flagDt)
    {
        numStates = hmm->numStates;
    }
    else
    {
        numStates = hmm->numStatesDt;
    }
    for(row = 0; row < numStates; row++)
    {
        for(col = 0; col < numStates-1; col++)
        {
            printf("%e,", hmm->pts[row][col]);
        }
        printf("%e\n", hmm->pts[row][numStates-1]);
    }
    return;
}

void Hmm_print_emissions(Hmm * hmm)
{
    assert(hmm);
    int i, j, ijidx;
    const int numStates = hmm->numStates;
    const int n = hmm->n;
    for(i = 0; i <= n; i++)
    {
        for(j = i; j <= n; j++)
        {
            ijidx = get_index(i, j, n);
            printf("%i\t%i\t%e\t%e\t%e\t%e\n", i, j, hmm->emissions[ijidx][0], 
                    hmm->emissions[ijidx][1], hmm->emissions[ijidx][2],
                    hmm->emissions[ijidx][3]);
        }
    }
    return;
}
