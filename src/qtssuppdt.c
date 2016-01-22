#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "definitions.h"
#include "hmm.h"
#include "qtscases.h"
#include "qtssuppdt.h"

double qts_case_A_supp(Hmm * hmm, int i, int j, int k, int l)
{
    assert(hmm);
    assert(0 <= i && i == k && k < j && j < l && l <= hmm->n);
    assert(hmm->flagDt);

    int a;
    const int n = hmm->n;
    const int rowIdx = get_index(i,j,n);
    double * const t = hmm->ts;
    double * const lam = hmm->lambdas;
    double * const del = hmm->deltas;
    double * const io = hmm->intervalOmegas;
    double ** const qts = hmm->qts;
    const double Es3 = hmm->Eijs3s[rowIdx];
    const double Es2 = hmm->Eijs2s[rowIdx];
    const double D3 = hmm->D3;
    const double lamd = hmm->lambdaDt;
    assert(D3 <= 0.0);

    double prob = 1.0/(2.0*Es2+Es3-D3) * lamd * (1.0 - exp(D3/lamd)) *
        exp(-3*get_omega_interval_Es3(hmm, 0, i, j)) * 
        exp(-2*get_omega_Es3_Es2(hmm, i, j)) * 
        exp(-get_omega_Es2_interval(hmm, l, i, j)) *
        ((l == n) ? 1.0 : (1.0 - exp(-io[l])));

    assert(0 <= prob && prob <= 1);
    return prob;
}

double qts_case_B_supp(Hmm * hmm, int i, int j, int k, int l)
{
    assert(hmm);
    assert(0 <= i && i == k && k < l && l < j && j <= hmm->n);
    assert(hmm->flagDt);

    int a;
    const int n = hmm->n;
    const int rowIdx = get_index(i,j,n);
    double * const t = hmm->ts;
    double * const lam = hmm->lambdas;
    double * const del = hmm->deltas;
    double * const io = hmm->intervalOmegas;
    double ** const qts = hmm->qts;
    const double Es3 = hmm->Eijs3s[rowIdx];
    const double Es2 = hmm->Eijs2s[rowIdx];
    const double D3 = hmm->D3;
    const double lamd = hmm->lambdaDt;

    double prob;

    prob = 1.0/(2.0*Es2+Es3-D3) * lamd * (1.0 - exp(D3/lamd)) * 
        exp(-3*get_omega_interval_Es3(hmm, 0, i, j)) * 
        exp(-2*get_omega_Es3_interval(hmm, l, i, j)) * 
        1.0/2.0 * (1.0 - exp(-2.0*io[l]));

    assert(0 <= prob && prob <= 1);
    return prob;
}

double qts_case_C_supp(Hmm * hmm, int i, int j, int k, int l)
{
    assert(hmm);
    // two cases for cases C and D
    assert( (0 <= k && k < i && i == l && l < j && j <= hmm->n) ||
            (0 <= k && k < i && i < j && j == l && l <= hmm->n) );
    assert(hmm->flagDt);

    int a;
    const int n = hmm->n;
    const int rowIdx = get_index(i,j,n);
    double * const t = hmm->ts;
    double * const lam = hmm->lambdas;
    double * const del = hmm->deltas;
    double * const io = hmm->intervalOmegas;
    double ** const qts = hmm->qts;
    const double Es3 = hmm->Eijs3s[rowIdx];
    const double Es2 = hmm->Eijs2s[rowIdx];
    const double D3 = hmm->D3;
    const double lamd = hmm->lambdaDt;

    double prob;

    prob = (lamd * (1.0 - exp(D3/lamd))) / (2.0*Es2+Es3-D3) * 2.0/3.0 *
        exp(-3.0*get_omega_interval_interval(hmm, 0, k)) *
        (1.0 - exp(-3.0*io[k]));

    assert(0 <= prob && prob <= 1);
    return prob;
}

double qts_case_D_supp(Hmm * hmm, int i, int j, int k, int l)
{
    assert(hmm);
    assert(0 <= k && k < i && i < j && j == l && l <= hmm->n);
    assert(hmm->flagDt);

    int a;
    const int n = hmm->n;
    const int rowIdx = get_index(i,j,n);
    double * const t = hmm->ts;
    double * const lam = hmm->lambdas;
    double * const del = hmm->deltas;
    double * const io = hmm->intervalOmegas;
    double ** const qts = hmm->qts;

    const double Es3 = hmm->Eijs3s[rowIdx];
    const double Es2 = hmm->Eijs2s[rowIdx];
    const double D3 = hmm->D3;
    const double lamd = hmm->lambdaDt;
    
    double prob;

    prob = (lamd * (1.0-exp(D3/lamd))) / (2.0*Es2+Es3-D3) * 2.0/3.0 *
        exp(-3.0*get_omega_interval_interval(hmm, 0, k)) *
        (1.0 - exp(-3.0*io[k]));

    assert(0 <= prob && prob <= 1);
    return prob;
}

double qts_case_E_supp(Hmm * hmm, int i, int j, int k, int l)
{
    assert(hmm);
    assert(0 <= i && i < k && k < j && j == l && l <= hmm->n);
    assert(hmm->flagDt);

    int a;
    const int n = hmm->n;
    const int rowIdx = get_index(i,j,n);
    double * const t = hmm->ts;
    double * const lam = hmm->lambdas;
    double * const del = hmm->deltas;
    double * const io = hmm->intervalOmegas;
    double ** const qts = hmm->qts;

    const double Es3 = hmm->Eijs3s[rowIdx];
    const double Es2 = hmm->Eijs2s[rowIdx];
    const double D3 = hmm->D3;
    const double lamd = hmm->lambdaDt;

    double prob;

    prob = (lamd * (1.0-exp(D3/lamd))) / (2.0*Es2+Es3-D3) *
        exp(-3.0*get_omega_interval_Es3(hmm, 0, i, j)) * 
        exp(-2.0*get_omega_Es3_interval(hmm, k, i, j)) *
        (1.0 - exp(-2.0*io[k]));

    assert(0 <= prob && prob <= 1);
    return prob;
}

double qts_case_F_supp(Hmm * hmm, int i, int j, int k, int l)
{
    assert(hmm);
    assert(0 <= i && i < k && j == k && k < l && l <= hmm->n);
    assert(hmm->flagDt);

    int a;
    const int n = hmm->n;
    const int rowIdx = get_index(i,j,n);
    double * const t = hmm->ts;
    double * const lam = hmm->lambdas;
    double * const del = hmm->deltas;
    double * const io = hmm->intervalOmegas;
    double ** const qts = hmm->qts;
    const double Es3 = hmm->Eijs3s[rowIdx];
    const double Es2 = hmm->Eijs2s[rowIdx];
    const double D3 = hmm->D3;
    const double lamd = hmm->lambdaDt;

    double prob;

    prob = (lamd * (1.0-exp(D3/lamd))) / (2.0*Es2+Es3-D3) *
        exp(-3.0*get_omega_interval_Es3(hmm, 0, i, j)) *
        exp(-2.0*get_omega_Es3_Es2(hmm, i, j)) *
        exp(-get_omega_Es2_interval(hmm, l, i, j)) *
        ((l == n) ? 1.0 :(1.0 - exp(-io[l])));

    assert(0 <= prob && prob <= 1);
    return prob;
}

double qts_case_G_supp(Hmm * hmm, int i, int j, int k, int l)
{
    assert(hmm);
    assert(0 <= i && i == k && k == l && l < j && j <= hmm->n);
    assert(hmm->flagDt);

    int a;
    const int n = hmm->n;
    const int rowIdx = get_index(i,j,n);
    double * const t = hmm->ts;
    double * const lam = hmm->lambdas;
    double * const del = hmm->deltas;
    double * const io = hmm->intervalOmegas;
    double ** const qts = hmm->qts;
    const double Es3 = hmm->Eijs3s[rowIdx];
    const double Es2 = hmm->Eijs2s[rowIdx];
    const double D3 = hmm->D3;
    const double lamd = hmm->lambdaDt;

    double prob;

    prob = (lamd * (1.0-exp(D3/lamd))) / (2.0*Es2+Es3-D3) *
        exp(-3.0*get_omega_interval_Es3(hmm, 0, i, j)) * 
        1.0/2.0 * (1.0 - exp(-2.0*(t[l+1] - Es3)/lam[l]));

    assert(0 <= prob && prob <= 1);
    return prob;
}

double qts_case_G2_supp(Hmm * hmm, int i, int j, int k, int l)
{
    assert(hmm);
    assert(0 <= i && i == k && k == l && l < j && j <= hmm->n);
    assert(hmm->flagDt);

    int a;
    const int n = hmm->n;
    const int rowIdx = get_index(i,j,n);
    double * const t = hmm->ts;
    double * const lam = hmm->lambdas;
    double * const del = hmm->deltas;
    double * const io = hmm->intervalOmegas;
    double ** const qts = hmm->qts;
    const double Es3 = hmm->Eijs3s[rowIdx];
    const double Es2 = hmm->Eijs2s[rowIdx];
    const double D3 = hmm->D3;
    const double lamd = hmm->lambdaDt;

    double prob;

    prob = (lamd * (1.0-exp(D3/lamd))) / (2.0*Es2+Es3-D3) * 2.0/3.0 *
        exp(-3.0*get_omega_interval_interval(hmm, 0, k)) *
        (1.0 - exp(-3*(Es3-t[k])/lam[k]));

    assert(0 <= prob && prob <= 1);
    return prob;
}

double qts_case_H_supp(Hmm * hmm, int i, int j, int k, int l)
{
    assert(hmm);
    assert(0 <= i && i < k && k == l && l == j && j <= hmm->n);
    assert(hmm->flagDt);

    int a;
    const int n = hmm->n;
    const int rowIdx = get_index(i,j,n);
    double * const t = hmm->ts;
    double * const lam = hmm->lambdas;
    double * const del = hmm->deltas;
    double * const io = hmm->intervalOmegas;
    double ** const qts = hmm->qts;
    const double Es3 = hmm->Eijs3s[rowIdx];
    const double Es2 = hmm->Eijs2s[rowIdx];
    const double D3 = hmm->D3;
    const double lamd = hmm->lambdaDt;

    double prob;

    prob = (lamd * (1.0-exp(D3/lamd))) / (2.0*Es2+Es3-D3) * 
        exp(-3.0*get_omega_interval_Es3(hmm, 0, i, j)) * 
        exp(-2.0*get_omega_Es3_interval(hmm, k, i, j)) *
        (1.0 - exp(-2.0*(Es2-t[k])/lam[k]));

    assert(0 <= prob && prob <= 1);
    return prob;
}

double qts_case_H2_supp(Hmm * hmm, int i, int j, int k, int l)
{
    assert(hmm);
    assert(0 <= i && i < k && k == l && l == j && j <= hmm->n);
    assert(hmm->flagDt);

    int a;
    const int n = hmm->n;
    const int rowIdx = get_index(i,j,n);
    double * const t = hmm->ts;
    double * const lam = hmm->lambdas;
    double * const del = hmm->deltas;
    double * const io = hmm->intervalOmegas;
    double ** const qts = hmm->qts;
    const double Es3 = hmm->Eijs3s[rowIdx];
    const double Es2 = hmm->Eijs2s[rowIdx];
    const double D3 = hmm->D3;
    const double lamd = hmm->lambdaDt;

    double prob;

    prob = (lamd * (1.0-exp(D3/lamd))) / (2.0*Es2+Es3-D3) *
        exp(-3.0*get_omega_interval_Es3(hmm, 0, i, j)) *
        exp(-2.0*get_omega_Es3_Es2(hmm, i, j)) *
        ((l == n) ? 1.0 : (1.0 - exp(-(t[l+1] - Es2)/lam[l])));

    assert(0 <= prob && prob <= 1);
    return prob;
}

double qts_case_I_supp(Hmm * hmm, int i, int j, int k, int l)
{
    assert(hmm);
    assert(0 <= i && i == j && j == k && k < l && l <= hmm->n);
    assert(hmm->flagDt);

    int a;
    const int n = hmm->n;
    const int rowIdx = get_index(i,j,n);
    double * const t = hmm->ts;
    double * const lam = hmm->lambdas;
    double * const del = hmm->deltas;
    double * const io = hmm->intervalOmegas;
    double ** const qts = hmm->qts;
    const double Es3 = hmm->Eijs3s[rowIdx];
    const double Es2 = hmm->Eijs2s[rowIdx];
    const double D3 = hmm->D3;
    const double lamd = hmm->lambdaDt;

    double prob;

    assert(l <= n);

    prob = (lamd * (1.0-exp(D3/lamd))) / (2.0*Es2+Es3-D3) * 
        exp(-3.0*get_omega_interval_Es3(hmm, 0, i, j)) *
        exp(-2.0*get_omega_Es3_Es2(hmm, i, j)) *
        exp(-get_omega_Es2_interval(hmm, l, i, j)) *
        ((l == n) ? 1.0 : (1.0 - exp(-io[l])));

    assert(0 <= prob && prob <= 1);
    return prob;
}

double qts_case_I2_supp(Hmm * hmm, int i, int j, int k, int l)
{
    assert(hmm);
    assert(0 <= i && i == j && j == k && k < l && l <= hmm->n);
    assert(hmm->flagDt);

    int a;
    const int n = hmm->n;
    const int rowIdx = get_index(i,j,n);
    double * const t = hmm->ts;
    double * const lam = hmm->lambdas;
    double * const del = hmm->deltas;
    double * const io = hmm->intervalOmegas;
    double ** const qts = hmm->qts;
    const double Es3 = hmm->Eijs3s[rowIdx];
    const double Es2 = hmm->Eijs2s[rowIdx];
    const double D3 = hmm->D3;
    const double lamd = hmm->lambdaDt;

    double prob;

    prob = (lamd * (1.0-exp(D3/lamd))) / (2.0*Es2+Es3-D3) *
        exp(-3.0*get_omega_interval_Es3(hmm, 0, i ,j)) *
        exp(-2.0*get_omega_Es3_Es2(hmm, i, j)) *
        exp(-get_omega_Es2_interval(hmm, l, i, j)) *
        ((l == n) ? 1.0 : (1.0 - exp(-io[l])));

    assert(0 <= prob && prob <= 1);
    return prob;
}

double qts_case_J_supp(Hmm * hmm, int i, int j, int k, int l)
{
    assert(hmm);
    assert(0 <= k && k < i && i == j && j == l && l <= hmm->n);
    assert(hmm->flagDt);

    int a;
    const int n = hmm->n;
    const int rowIdx = get_index(i,j,n);
    double * const t = hmm->ts;
    double * const lam = hmm->lambdas;
    double * const del = hmm->deltas;
    double * const io = hmm->intervalOmegas;
    double ** const qts = hmm->qts;
    const double Es3 = hmm->Eijs3s[rowIdx];
    const double Es2 = hmm->Eijs2s[rowIdx];
    const double D3 = hmm->D3;
    const double lamd = hmm->lambdaDt;

    double prob;

    prob = (lamd * (1.0-exp(D3/lamd))) / (2.0*Es2+Es3-D3) * 2.0/3.0 *
        exp(-3.0*get_omega_interval_interval(hmm, 0, k)) *
        (1.0 - exp(-3.0*io[k]));

    assert(0 <= prob && prob <= 1);
    return prob;
}

double qts_case_J2_supp(Hmm * hmm, int i, int j, int k, int l)
{
    assert(hmm);
    assert(0 <= k && k < i && i == j && j == l && l <= hmm->n);
    assert(hmm->flagDt);

    int a;
    const int n = hmm->n;
    const int rowIdx = get_index(i,j,n);
    double * const t = hmm->ts;
    double * const lam = hmm->lambdas;
    double * const del = hmm->deltas;
    double * const io = hmm->intervalOmegas;
    double ** const qts = hmm->qts;
    const double Es3 = hmm->Eijs3s[rowIdx];
    const double Es2 = hmm->Eijs2s[rowIdx];
    const double D3 = hmm->D3;
    const double lamd = hmm->lambdaDt;

    double prob;

    prob = (lamd * (1.0-exp(D3/lamd))) / (2.0*Es2+Es3-D3) * 2.0/3.0 *
        exp(-3.0*get_omega_interval_interval(hmm, 0, k)) *
        (1.0 - exp(-3.0*io[k]));

    assert(0 <= prob && prob <= 1);
    return prob;
}
