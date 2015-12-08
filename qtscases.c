#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "hmm.h"
#include "qtscases.h"

double qts_case_A(Hmm * hmm, int i, int j, int k, int l)
{
    assert(hmm);
    assert(0 <= i && i == k && k < j && j < l && l <= hmm->n);

    int a;
    const int n = hmm->n;
    const int rowIdx = get_index(i,j,n);
    const int colIdx = get_index(k,l,n);
    double * const t = hmm->ts;
    double * const lam = hmm->lambdas;
    double * const del = hmm->deltas;
    double * const io = hmm->intervalOmegas;
    double * const pi = hmm->pis;
    double * const es3 = hmm->Eijs3s;
    double * const es2 = hmm->Eijs2s;
    double ** const qts = hmm->qts;

    double prob;

    double sum1 = 0.0;
    for(a = 0; a <= i-1; a++) // note <= here and below
    {
        sum1 += exp(-3*get_omega_interval_Es3(hmm, a+1, i, j)) *
            (1-exp(-3*io[a])) * lam[a]/3.0;
    }

    double sum2 = 0.0;
    for(a = i+1; a <= j-1; a++)
    {
        sum2 += exp(-2*get_omega_interval_Es2(hmm, a+1, i, j)) *
            (1.0 - exp(-2*io[a])) * lam[a]/2.0;
    }

    prob = 1.0/(2*es2[rowIdx]+es3[rowIdx])*exp(-2*get_omega_Es3_Es2(hmm,i,j)) /
        lam[l] * (sum1 + (1-exp(-3*(es3[rowIdx]-t[i])/lam[i])) * lam[i]/3.0) * 
        exp(-get_omega_Es2_interval(hmm, l, i, j)) * lam[l] *
        ((l == n) ? 1.0 : (1.0 - exp(-io[l]))) + 
        2.0/(2.0*es2[rowIdx]+es3[rowIdx]) / lam[l] *
        ( exp(-2*get_omega_interval_Es2(hmm, i+1, i, j)) *
           (1.0-exp(-2*(t[i+1]-es3[rowIdx])/lam[i])) * lam[i]/2.0 + sum2 +
           (1.0-exp(-2*(es2[rowIdx]-t[j])/lam[j])) * lam[j]/2.0 ) *
        exp(-get_omega_Es2_interval(hmm, l, i, j)) * lam[l] *
        ((l == n) ? 1.0 : (1.0 - exp(-io[l]))); // note the 'delta' function here

    assert(0 <= prob && prob <= 1);
    return prob;
}

double qts_case_B(Hmm * hmm, int i, int j, int k, int l)
{
    assert(hmm);
    assert(0 <= i && i == k && k < l && l < j && j <= hmm->n);

    int a;
    const int n = hmm->n;
    const int rowIdx = get_index(i,j,n);
    const int colIdx = get_index(k,l,n);
    double * const t = hmm->ts;
    double * const lam = hmm->lambdas;
    double * const del = hmm->deltas;
    double * const io = hmm->intervalOmegas;
    double * const pi = hmm->pis;
    double * const es3 = hmm->Eijs3s;
    double * const es2 = hmm->Eijs2s;
    double ** const qts = hmm->qts;

    double prob;

    double sum1 = 0.0;
    for(a = 0; a <= i-1; a++)
    {
        sum1 += exp(-3*get_omega_interval_Es3(hmm, a+1, i, j)) *
            (1.0 - exp(-3*io[a])) * lam[a]/3.0;
    }

    double sum2 = 0.0;
    for(a = i+1; a <= l-1; a++)
    {
        sum2 += (1-exp(-2*io[a]))* lam[a]/2.0 * 
            exp(-2*get_omega_interval_interval(hmm, a+1, l));
        // note that we are not including the lam[l]/2.0 * (1-...)
        // factor, as it doesn't depend on a. will be included later.
    }

    prob = 1.0/(2.0*es2[rowIdx]+es3[rowIdx])*1.0/lam[l] * 
        (sum1 + (1.0 - exp(-3*(es3[rowIdx]-t[i])/lam[i])) * lam[i]/3.0) *
        exp(-get_omega_Es3_interval(hmm, l, i, j)) * lam[l]/2.0 *
        (1.0-exp(-2*io[l])) +
        2.0/(2*es2[rowIdx]+es3[rowIdx])*1.0/lam[l] *
        (sum2 * lam[l]/2.0 * (1-exp(-2*io[l]))) +
        2.0/(2*es2[rowIdx]+es3[rowIdx])*1.0/lam[l] * lam[l]/2.0 *
        (del[l] - lam[l]/2.0 * (1-exp(-2*io[l])));

    assert(0 <= prob && prob <= 1);
    return prob;
}

double qts_case_C(Hmm * hmm, int i, int j, int k, int l)
{
    assert(hmm);
    // two cases for cases C and D
    assert( (0 <= k && k < i && i == l && l < j && j <= hmm->n) ||
            (0 <= k && k < i && i < j && j == l) );

    int a;
    const int n = hmm->n;
    const int rowIdx = get_index(i,j,n);
    const int colIdx = get_index(k,l,n);
    double * const t = hmm->ts;
    double * const lam = hmm->lambdas;
    double * const del = hmm->deltas;
    double * const io = hmm->intervalOmegas;
    double * const pi = hmm->pis;
    double * const es3 = hmm->Eijs3s;
    double * const es2 = hmm->Eijs2s;
    double ** const qts = hmm->qts;

    double prob;

    double sum1 = 0.0;
    for(a = 0; a <= k-1; a++)
    {
        sum1 += (1.0 - exp(-3*io[a])) * lam[a]/3.0 *
            exp(-3*get_omega_interval_interval(hmm, a+1, k));
    }

    prob = 1.0/(2.0*es2[rowIdx]+es3[rowIdx])* 2.0/lam[k] * 
        (lam[k]/3.0 * (1-exp(-3*io[k])) * sum1 + 
         lam[k]/3.0 * (del[k] - lam[k]/3.0 * (1-exp(-3*io[k]))));

    assert(0 <= prob && prob <= 1);
    return prob;
}

double qts_case_D(Hmm * hmm, int i, int j, int k, int l)
{
    assert(hmm);
    assert(0 <= k && k < i && i < j && j == l && l <= hmm->n);

    // case D is 2x case C
    double prob = 2.0*qts_case_C(hmm, i, j, k, l);

    assert(0 <= prob && prob <= 1);
    return prob;
}

double qts_case_E(Hmm * hmm, int i, int j, int k, int l)
{
    assert(hmm);
    assert(0 <= i && i < k && k < j && j == l && l <= hmm->n);

    int a;
    const int n = hmm->n;
    const int rowIdx = get_index(i,j,n);
    const int colIdx = get_index(k,l,n);
    double * const t = hmm->ts;
    double * const lam = hmm->lambdas;
    double * const del = hmm->deltas;
    double * const io = hmm->intervalOmegas;
    double * const pi = hmm->pis;
    double * const es3 = hmm->Eijs3s;
    double * const es2 = hmm->Eijs2s;
    double ** const qts = hmm->qts;

    double prob;

    double sum1 = 0.0;
    for(a = 0; a <= i-1; a++)
    {
        sum1 += exp(-3*get_omega_interval_Es3(hmm, a+1, i, j)) *
            (1.0-exp(-3*io[a])) * lam[a]/3.0;
    }

    prob = 2.0/(2.0*es2[rowIdx]+es3[rowIdx]) * 2.0/lam[k] * 
        (sum1 + (1-exp(-3*(es3[rowIdx]-t[i])/lam[i])) * lam[i]/3.0) *
        exp(-2*get_omega_Es3_interval(hmm, k, i, j)) * lam[k]/2.0 *
        (1.0 - exp(-2*io[k]));

    assert(0 <= prob && prob <= 1);
    return prob;
}

double qts_case_F(Hmm * hmm, int i, int j, int k, int l)
{
    assert(hmm);
    assert(0 <= i && i < k && j == k && k < l && l <= hmm->n);

    int a;
    const int n = hmm->n;
    const int rowIdx = get_index(i,j,n);
    const int colIdx = get_index(k,l,n);
    double * const t = hmm->ts;
    double * const lam = hmm->lambdas;
    double * const del = hmm->deltas;
    double * const io = hmm->intervalOmegas;
    double * const pi = hmm->pis;
    double * const es3 = hmm->Eijs3s;
    double * const es2 = hmm->Eijs2s;
    double ** const qts = hmm->qts;

    double prob;

    double sum1 = 0.0;
    for(a = 0; a <= i-1; a++)
    {
        sum1 += exp(-3*get_omega_interval_Es3(hmm, a+1, i, j)) * 
            (1.0 - exp(-3*io[a]));
    }

    prob = 2.0/(2.0*es2[rowIdx]+es3[rowIdx]) *
        exp(-2*get_omega_Es3_Es2(hmm, i, j)) / lam[l] *
        (sum1 + (1.0 - exp(-3*(es3[rowIdx] - t[i])/lam[i])) * lam[i]/3.0) *
        exp(-get_omega_Es2_interval(hmm, l, i, j)) * lam[l] *
        ((l == n) ? 1.0 : 1.0-exp(-io[l])); // 'delta' func again w/ cond. exp.

    assert(0 <= prob && prob <= 1);
    return prob;
}

double qts_case_G(Hmm * hmm, int i, int j, int k, int l)
{
    assert(hmm);
    assert(0 <= i && i == k && k == l && l < j && j <= hmm->n);

    int a;
    const int n = hmm->n;
    const int rowIdx = get_index(i,j,n);
    const int colIdx = get_index(k,l,n);
    double * const t = hmm->ts;
    double * const lam = hmm->lambdas;
    double * const del = hmm->deltas;
    double * const io = hmm->intervalOmegas;
    double * const pi = hmm->pis;
    double * const es3 = hmm->Eijs3s;
    double * const es2 = hmm->Eijs2s;
    double ** const qts = hmm->qts;

    double prob;

    double sum1 = 0.0;
    for(a = 0; a <= i-1; a++)
    {
        sum1 += exp(-3*get_omega_interval_Es3(hmm, a+1, i, j)) *
            (1.0 - exp(-3*io[a])) * lam[a]/3.0;
    }

    prob = 1.0/(2.0*es3[rowIdx]+es3[rowIdx]) * ( 1.0/lam[k] *
            (sum1 + (1-exp(-3*(es3[rowIdx]-t[i])/lam[i])) * lam[i]/3.0) *
            (1-exp(-2*(t[k+1]-es3[rowIdx])/lam[k])) +
            t[k+1] - es3[rowIdx] - 
            lam[k]/2.0 *(1-exp(-2*(t[k+1]-es3[rowIdx])/lam[k])));

    assert(0 <= prob && prob <= 1);
    return prob;
}

double qts_case_G2(Hmm * hmm, int i, int j, int k, int l)
{
    assert(hmm);
    assert(0 <= i && i == k && k == l && l < j && j <= hmm->n);

    int a;
    const int n = hmm->n;
    const int rowIdx = get_index(i,j,n);
    const int colIdx = get_index(k,l,n);
    double * const t = hmm->ts;
    double * const lam = hmm->lambdas;
    double * const del = hmm->deltas;
    double * const io = hmm->intervalOmegas;
    double * const pi = hmm->pis;
    double * const es3 = hmm->Eijs3s;
    double * const es2 = hmm->Eijs2s;
    double ** const qts = hmm->qts;

    double prob;

    double sum1 = 0.0;
    for(a = 0; a <= k-1; a++)
    {
        sum1 += (1-exp(-3*io[a]))*lam[a]/3.0*
            exp(-3*get_omega_interval_interval(hmm, a+1, k));
    }

    prob = 1.0/(2.0*es2[rowIdx] + es3[rowIdx])*2.0/lam[k] * 
        (sum1 * lam[k]/3.0 *
        (1.0-exp(-3*(es3[rowIdx]-t[k])/lam[k])) + 
        lam[k]/3.0 * (es3[rowIdx]-t[k]-lam[k]/3.0 * 
            (1-exp(-3*(es3[rowIdx]-t[k])/lam[k]))));

    assert(0 <= prob && prob <= 1);
    return prob;
}

double qts_case_H(Hmm * hmm, int i, int j, int k, int l)
{
    assert(hmm);
    assert(0 <= i && i < k && k == l && l == j && j <= hmm->n);

    int a;
    const int n = hmm->n;
    const int rowIdx = get_index(i,j,n);
    const int colIdx = get_index(k,l,n);
    double * const t = hmm->ts;
    double * const lam = hmm->lambdas;
    double * const del = hmm->deltas;
    double * const io = hmm->intervalOmegas;
    double * const pi = hmm->pis;
    double * const es3 = hmm->Eijs3s;
    double * const es2 = hmm->Eijs2s;
    double ** const qts = hmm->qts;

    double prob;

    double sum1 = 0.0;
    for(a = 0; a <= i-1; a++)
    {
        sum1 += exp(-3*get_omega_interval_Es3(hmm, a+1, i, j))
            * (1.0 - exp(-3*io[a])) * lam[a]/3.0;
    }

    prob = 1.0/(2.0*es2[rowIdx]+es3[rowIdx])*4.0/lam[k] *
        (sum1 + (1-exp(-3*(es3[rowIdx]-t[i])/lam[i]))*lam[i]/3) *
        lam[k]/2.0*(1-exp(-2*(es2[rowIdx]-t[k])/lam[k]));

    assert(0 <= prob && prob <= 1);
    return prob;
}

double qts_case_H2(Hmm * hmm, int i, int j, int k, int l)
{
    assert(hmm);
    assert(0 <= i && i < k && k == l && l == j && j <= hmm->n);

    int a;
    const int n = hmm->n;
    const int rowIdx = get_index(i,j,n);
    const int colIdx = get_index(k,l,n);
    double * const t = hmm->ts;
    double * const lam = hmm->lambdas;
    double * const del = hmm->deltas;
    double * const io = hmm->intervalOmegas;
    double * const pi = hmm->pis;
    double * const es3 = hmm->Eijs3s;
    double * const es2 = hmm->Eijs2s;
    double ** const qts = hmm->qts;

    double prob;

    double sum1 = 0.0;
    for(a = 0; a <= i-1; a++)
    {
        sum1 += exp(-3*get_omega_interval_Es3(hmm, a+1, i, j))
            * (1.0 - exp(-3*io[a])) * lam[a]/3.0;
    }

    prob = 2.0/(2.0*es2[rowIdx]+es3[rowIdx]) / lam[k] * 
        exp(-2*get_omega_Es3_Es2(hmm, i, j)) *
        (sum1 + (1-exp(-3*(es3[rowIdx]-t[i])/lam[i])) * lam[i]/3.0) *
        lam[k] * ((k == n) ? 1.0 : 1.0-exp(-(t[k+1]-es2[rowIdx])/lam[k]));

    assert(0 <= prob && prob <= 1);
    return prob;
}

double qts_case_I(Hmm * hmm, int i, int j, int k, int l)
{
    assert(hmm);
    assert(0 <= i && i == j && j == k && k < l && l <= hmm->n);

    int a;
    const int n = hmm->n;
    const int rowIdx = get_index(i,j,n);
    const int colIdx = get_index(k,l,n);
    double * const t = hmm->ts;
    double * const lam = hmm->lambdas;
    double * const del = hmm->deltas;
    double * const io = hmm->intervalOmegas;
    double * const pi = hmm->pis;
    double * const es3 = hmm->Eijs3s;
    double * const es2 = hmm->Eijs2s;
    double ** const qts = hmm->qts;

    double prob;

    double sum1 = 0.0;
    for(a = 0; a <= i-1; a++)
    {
        sum1 += exp(-3*get_omega_interval_Es3(hmm, a+1, i, j))
            * (1.0 - exp(-3*io[a])) * lam[a]/3.0;
    }

    // TODO: double check conditional expressions
    prob = 1.0/(2.0*es2[rowIdx]+es3[rowIdx]) * 
        (exp(-2*get_omega_Es3_Es2(hmm, i, j)) / lam[l] *
         ( sum1 + (1-exp(-3*(es3[rowIdx]-t[i])/lam[i]))*lam[i]/3.0) * 
         exp(-get_omega_Es2_interval(hmm, l, i, j)) * lam[l] *
         ((l == n) ? 1 : 1-exp(-io[l])) +
        2.0/lam[l] * lam[i]/2.0 *
        (1-exp(-2*(es2[rowIdx]-es3[rowIdx])/lam[i])) * 
        exp(-get_omega_Es2_interval(hmm, l, i, j)) * lam[l] * 
        ((l == n) ? 1 : 1-exp(-io[l])));

    assert(0 <= prob && prob <= 1);
    return prob;
}

double qts_case_I2(Hmm * hmm, int i, int j, int k, int l)
{
    assert(hmm);
    assert(0 <= i && i == j && j == k && k < l && l <= hmm->n);

    int a;
    const int n = hmm->n;
    const int rowIdx = get_index(i,j,n);
    const int colIdx = get_index(k,l,n);
    double * const t = hmm->ts;
    double * const lam = hmm->lambdas;
    double * const del = hmm->deltas;
    double * const io = hmm->intervalOmegas;
    double * const pi = hmm->pis;
    double * const es3 = hmm->Eijs3s;
    double * const es2 = hmm->Eijs2s;
    double ** const qts = hmm->qts;

    double prob;

    double sum1 = 0.0;
    for(a = 0; a <= i-1; a++)
    {
        sum1 += exp(-3*get_omega_interval_Es3(hmm, a+1, i, j))
            * (1.0 - exp(-3*io[a])) * lam[a]/3.0;
    }

    // TODO: double check conditional expressions
    prob = 2.0/(2.0*es2[rowIdx]+es3[rowIdx]) * 
        exp(-2*get_omega_Es3_Es2(hmm, i, j)) / lam[l] *
         exp(-get_omega_Es2_interval(hmm, l, i, j)) * lam[l] *
         ((l == n) ? 1 : 1-exp(-io[l])) * 
         ( sum1 + (1-exp(-3*(es3[rowIdx]-t[i])/lam[i]))*lam[i]/3.0);

    assert(0 <= prob && prob <= 1);
    return prob;
}

double qts_case_J(Hmm * hmm, int i, int j, int k, int l)
{
    assert(hmm);
    assert(0 <= k && k < i && i == j && j == l && l <= hmm->n);

    int a;
    const int n = hmm->n;
    const int rowIdx = get_index(i,j,n);
    const int colIdx = get_index(k,l,n);
    double * const t = hmm->ts;
    double * const lam = hmm->lambdas;
    double * const del = hmm->deltas;
    double * const io = hmm->intervalOmegas;
    double * const pi = hmm->pis;
    double * const es3 = hmm->Eijs3s;
    double * const es2 = hmm->Eijs2s;
    double ** const qts = hmm->qts;

    double prob;

    double sum1 = 0.0;
    for(a = 0; a <= k-1; a++)
    {
        sum1 += (1-exp(-3*io[a]))*lam[a]/3.0 * 
            exp(-3*get_omega_interval_interval(hmm, a+1, k));
    }

    // TODO: double check conditional expressions
    prob = 2.0/(2.0*es2[rowIdx]+es3[rowIdx]) * 2.0 / lam[k] *
        (lam[k]/3.0*(1-exp(-3*io[k])) * sum1 +
         lam[k]/3.0*(del[k]-lam[k]/3.0*(1-exp(-3*io[k]))));

    assert(0 <= prob && prob <= 1);
    return prob;
}

double qts_case_J2(Hmm * hmm, int i, int j, int k, int l)
{
    assert(hmm);
    assert(0 <= k && k < i && i == j && j == l && l <= hmm->n);

    int a;
    const int n = hmm->n;
    const int rowIdx = get_index(i,j,n);
    const int colIdx = get_index(k,l,n);
    double * const t = hmm->ts;
    double * const lam = hmm->lambdas;
    double * const del = hmm->deltas;
    double * const io = hmm->intervalOmegas;
    double * const pi = hmm->pis;
    double * const es3 = hmm->Eijs3s;
    double * const es2 = hmm->Eijs2s;
    double ** const qts = hmm->qts;

    double prob;

    double sum1 = 0.0;
    for(a = 0; a <= k-1; a++)
    {
        sum1 += (1-exp(-3*io[a]))*lam[a]/3.0 * 
            exp(-3*get_omega_interval_interval(hmm, a+1, k));
    }

    prob = 2.0/(2.0*es2[rowIdx]+es3[rowIdx]) * 1.0 / lam[k] *
        (lam[k]/3.0*(1-exp(-3*io[k])) * sum1 +
         lam[k]/3.0*(del[k]-lam[k]/3.0*(1-exp(-3*io[k]))));

    assert(0 <= prob && prob <= 1);
    return prob;
}
