#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "definitions.h"
#include "hmm.h"
#include "data.h"
#include "options.h"
#include "em.h"
#include "random.h"
#include "asa047.h"

Em em;

void Em_init(Em * em, Data * dat, double * ts,
        double initRho, int numFreeLambdas, 
        int n, int numEmIterations, int * lambdaCounts, int asexEnabled)
{
    assert(em && dat && ts);
    assert(dat->numSeqs > 0);
    em->dat = dat;
    em->numSeqs = dat->numSeqs;
    em->seqtype = dat->seqtype;
    em->asexEnabled = asexEnabled;
    assert(em->asexEnabled == 0 || em->asexEnabled == 1);
    em->lambdaCounts = lambdaCounts;
    em->forward = (double ***)chmalloc(sizeof(double **) * em->numSeqs);
    em->backward = (double ***)chmalloc(sizeof(double **) * em->numSeqs);
    em->normConst = (double **)chmalloc(sizeof(double *) * em->numSeqs);
    em->gamma = (double ***)chmalloc(sizeof(double **) * em->numSeqs);
    em->freeLambdas = (double *)chmalloc(sizeof(double) * numFreeLambdas);
    em->numFreeLambdas = numFreeLambdas;

    double thetaAndTd[2];
    Em_get_initial_theta_and_Td(em, &(thetaAndTd[0]));

    const double initTheta = thetaAndTd[0];
    const double initTd = thetaAndTd[1];
    DEBUGREPORTF(initTheta);
    DEBUGREPORTF(initTd);

    int i, j;

    for(i = 0; i < em->numFreeLambdas; i++)
    {
        em->freeLambdas[i] = 1.0;
    }

    double * lambdas = (double *)chmalloc(sizeof(double) * (n+1));
    for(i = 0; i < n+1; i++)
    {
        lambdas[i] = 1.0;
    }

    Hmm_init(&(em->hmm[0]), n, ts);
    Hmm_init(&(em->hmm[1]), n, ts);
    Hmm_make_hmm(&(em->hmm[0]), lambdas, n, initTheta, initRho, initTd);
    em->hmmFlag = 0;
    em->curIteration = 0;
    em->maxIterations = numEmIterations;

    const int numHmmStates = (n+1)*(n+2)/2;
    em->numHmmStates = numHmmStates;

    for(i = 0; i < em->numSeqs; i++)
    {
        em->forward[i] = (double **)chmalloc(sizeof(double *) * 
                dat->seqs[i].len);
        em->backward[i] = (double **)chmalloc(sizeof(double *) * 
                dat->seqs[i].len);
        em->gamma[i] = (double **)chmalloc(sizeof(double *) * 
                dat->seqs[i].len);
        em->normConst[i] = (double *)chmalloc(sizeof(double) * 
                dat->seqs[i].len);
        for(j = 0; j < dat->seqs[i].len; j++)
        {
            em->forward[i][j] = (double *)chmalloc(sizeof(double) * 
                    numHmmStates);
            em->backward[i][j] = (double *)chmalloc(sizeof(double) * 
                    numHmmStates);
            em->gamma[i][j] = (double *)chmalloc(sizeof(double) * 
                    numHmmStates);
        }
    }
    em->expectTransitions = (double **)chmalloc(sizeof(double *) * numHmmStates);
    for(i = 0; i < numHmmStates; i++)
    {
        em->expectTransitions[i] = (double *)chmalloc(sizeof(double) * numHmmStates);
    }
    em->expectEmissions = (fourd *)chmalloc(sizeof(fourd) * numHmmStates);
    free(lambdas);
    return;
}

void Em_free(Em * em)
{
    int i, j, k;
    for(i = 0; i < em->numSeqs; i++)
    {
        for(j = 0; j < em->dat->seqs[i].len; j++)
        {
            free(em->forward[i][j]);
            free(em->backward[i][j]);
            free(em->gamma[i][j]);
        }
        free(em->forward[i]);
        free(em->backward[i]);
        free(em->gamma[i]);
        free(em->normConst[i]);
    }
    for(i = 0; i < em->numHmmStates; i++)
    {
        free(em->expectTransitions[i]);
    }
    free(em->normConst);
    free(em->expectTransitions);
    free(em->forward);
    free(em->backward);
    free(em->gamma);
    free(em->freeLambdas);
    free(em->expectEmissions);
    Hmm_free(&(em->hmm[0]));
    Hmm_free(&(em->hmm[1]));
    return;
}

double Em_get_initial_rho(Data * dat)
{
    const int n = 8;
    double * ts = chmalloc(sizeof(double) * (n+1));
    get_ts_msmc(ts, n);
    int lambdaCounts = 8;
    Em tempEm;
    double initRhos[6] = {5e-3, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1};
    double maxloglike = -1.0*INFINITY;
    double loglike;
    int maxLLidx = 0;
    const int numInitRhos = 6;
    int i;
    for(i = 0; i < numInitRhos; i++)
    {
        Em_init(&tempEm, dat, ts, initRhos[i], 0, n, 2, &lambdaCounts, 1);
        Em_get_forward(&tempEm);
        Em_get_backward(&tempEm);
        Em_get_expectations(&tempEm);
        loglike = Em_get_loglikelihood(&tempEm);
        if(loglike > maxloglike)
        {
            maxLLidx = i;
            maxloglike = loglike;
        }
        Em_free(&tempEm);
    }
    free(ts);
    double initRho = initRhos[maxLLidx];
    return initRho;
}

void Em_get_initial_theta_and_Td(Em * em, double * out)
{
    assert(em && em->dat);
    int i, j, seqLen, totalLen = 0;
    int sfs[] = {0,0,0};
    char * seqData;
    for(i = 0; i < em->numSeqs; i++)
    {
        seqLen = em->dat->seqs[i].len;
        seqData = em->dat->seqs[i].data;
        for(j = 0; j < seqLen; j++)
        {
            totalLen++;
            if(seqData[j] > 0)
            {
                sfs[seqData[j]-1] += 1;
            }
        }
    }
    double thetaDoubleton = 2.0*(sfs[1]+sfs[2])/(double)totalLen;
    double thetaSingleton = (sfs[0]+sfs[2])/(double)totalLen;
    double meanTheta = (thetaDoubleton + thetaSingleton)/2.0;
    out[0] = meanTheta;

    double R = (double)(sfs[0]+sfs[2]) / (sfs[0] + sfs[1] + 2*sfs[2]);
    double initTd = (3.0*R-2.0) / (3.0*(1.0-R));
    if(initTd < 0)
    {
        initTd = 0.0;
    }
    out[1] = initTd;
    return;
}

void Em_get_forward(Em * em)
{
    const int hmmIdx = em->hmmFlag;
    double * const pis = em->hmm[hmmIdx].pis;
    double ** const pts = em->hmm[hmmIdx].pts;
    double *** const forward = em->forward;
    const int numEmissionStates = (em->seqtype == polarized) ? 4 : 2;
    const int numHmmStates = em->hmm[hmmIdx].numStates;
    assert(em->hmm[hmmIdx].numStates == em->hmm[!hmmIdx].numStates);
    const int numSeqs = em->numSeqs;
    fourd * const emissions = em->hmm[hmmIdx].emissions;

    int i, j, k, l;
    int seqLen;
    char * seqData;
    double ** seqFor;

    double sum1, sum2;

    for(i = 0; i < numSeqs; i++)
    {
        seqData = em->dat->seqs[i].data;
        seqLen = em->dat->seqs[i].len;
        seqFor = forward[i];
        // deal first with first position in sequence
        sum1 = 0.0;
        for(k = 0; k < numHmmStates; k++)
        {
            seqFor[0][k] = pis[k] * emissions[k][seqData[0]];
            assert(seqFor[0][k] > 0.0);
            sum1 += seqFor[0][k];
        }
        em->normConst[i][0] = sum1;
        for(k = 0; k < numHmmStates; k++)
        {
            seqFor[0][k] /= em->normConst[i][0];
        }

        for(j = 1; j < seqLen; j++)
        {
            assert(seqData[j] < numEmissionStates);
            sum2 = 0.0;
            for(k = 0; k < numHmmStates; k++)
            {
                sum1 = 0.0;
                for(l = 0; l < numHmmStates; l++)
                {
                    if(pts[l][k])
                    {
                        sum1 += seqFor[j-1][l] * pts[l][k];
                    }
                }
                seqFor[j][k] = emissions[k][seqData[j]] * sum1;
                sum2 += seqFor[j][k];
            }
            em->normConst[i][j] = sum2;
            for(k = 0; k < numHmmStates; k++)
            {
                seqFor[j][k] /= em->normConst[i][j];
            }
        }
    }
    return;
}

void Em_get_backward(Em * em)
{
    const int hmmIdx = em->hmmFlag;
    double * const pis = em->hmm[hmmIdx].pis;
    double ** const pts = em->hmm[hmmIdx].pts;
    double *** const backward = em->backward;
    const int numEmissionStates = (em->seqtype == polarized) ? 4 : 2;
    const int numHmmStates = em->hmm[hmmIdx].numStates;
    const int numSeqs = em->numSeqs;
    fourd * const emissions = em->hmm[hmmIdx].emissions;

    int i, j, k, l;
    int seqLen;
    char * seqData;
    double ** seqBack;

    double sum;

    for(i = 0; i < numSeqs; i++)
    {
        seqData = em->dat->seqs[i].data;
        seqLen = em->dat->seqs[i].len;
        seqBack = backward[i];
        for(k = 0; k < numHmmStates; k++)
        {
            seqBack[seqLen-1][k] = 1.0;
        }
        for(j = seqLen-2; j >= 0; j--)
        {
            assert(seqData[j] < numEmissionStates);
            for(k = 0; k < numHmmStates; k++)
            {
                sum = 0.0;
                for(l = 0; l < numHmmStates; l++)
                {
                    if(pts[k][l])
                    {
                        sum += seqBack[j+1][l] * pts[k][l] * emissions[l][seqData[j+1]];
                    }
                }
                assert(sum > 0);
                seqBack[j][k] = sum / em->normConst[i][j+1];
            }
        }
    }
    return;
}

void Em_get_expectations(Em * em)
{
    const int hmmIdx = em->hmmFlag;
    double *** const forward = em->forward;
    double *** const backward = em->backward;
    double *** const gamma = em->gamma;
    fourd * const expectEmissions = em->expectEmissions;
    double ** const expectTransitions = em->expectTransitions;
    double ** const pts = em->hmm[hmmIdx].pts;
    fourd * const emissions = em->hmm[hmmIdx].emissions;
    const int numHmmStates = em->hmm[hmmIdx].numStates;
    const int numSeqs = em->numSeqs;

    char * seqData;

    double ** seqFor, ** seqBack;

    int i, j, k, l, seqLen;

    double sum, thisExpect;
    
    // zero expectations
    for(j = 0; j < numHmmStates; j++)
    {
        for(k = 0; k < numHmmStates; k++)
        {
            expectTransitions[j][k] = 0.0;
        }
        for(k = 0; k < 4; k++)
        {
            expectEmissions[j][k] = 0.0;
        }
    }

    for(i = 0; i < numSeqs; i++)
    {
        seqData = em->dat->seqs[i].data;
        seqLen = em->dat->seqs[i].len;
        assert(seqLen > 0);

        seqFor = forward[i];
        seqBack = backward[i];

        // calculate gammas
        for(j = 0; j < seqLen; j++)
        {
            sum = 0.0;
            for(l = 0; l < numHmmStates; l++)
            {
                em->gamma[i][j][l] = seqFor[j][l]*seqBack[j][l];
                sum += em->gamma[i][j][l];
            }
            for(l = 0; l < numHmmStates; l++)
            {
                em->gamma[i][j][l] /= sum;
                expectEmissions[l][seqData[j]] += em->gamma[i][j][l];
                assert(0 < em->gamma[i][j][l] && em->gamma[i][j][l] < 1);
            }
        }

        // now calculate expected number of transitions
        for(j = 0; j < seqLen-1; j++)
        {
            for(k = 0; k < numHmmStates; k++)
            {
                for(l = 0; l < numHmmStates; l++)
                {
                    if(pts[k][l])
                    {
                        thisExpect = em->gamma[i][j][k] * pts[k][l] * 
                            seqBack[j+1][l] * emissions[l][seqData[j+1]] / 
                            (em->normConst[i][j+1] * seqBack[j][k]);
                        expectTransitions[k][l] += thisExpect;
                        assert(0 <= thisExpect && thisExpect <= 1);
                    }
                }
            }
        }
    }
    return;
}

double Em_get_expected_log_likelihood(Em * em, const int hmmIdx)
{
    Hmm * const hmm = &(em->hmm[hmmIdx]);
    double *** const gamma = em->gamma;
    double ** const expectTransitions = em->expectTransitions;
    fourd * expectEmissions = em->expectEmissions;
    double ** const pts = hmm->pts;
    double * const pis = hmm->pis;
    fourd * const emissions = hmm->emissions;
    const int numHmmStates = hmm->numStates;
    const int numSeqs = em->numSeqs;

    int i,j,k,l, seqLen;

    char * seqData;

    double loglike = 0.0;
    for(i = 0; i < numSeqs; i++)
    {
        seqData = em->dat->seqs[i].data;
        seqLen = em->dat->seqs[i].len;
        for(j = 0; j < numHmmStates; j++)
        {
            loglike += gamma[i][0][j] * log(pis[j]);
        }
    }
    for(i = 0; i < numHmmStates; i++)
    {
        for(j = 0; j < 4; j++)
        {
            loglike += log(emissions[i][j]) * expectEmissions[i][j];
        }
        for(j = 0; j < numHmmStates; j++)
        {
            if(pts[i][j] || expectTransitions[i][j] > 0.0)
            {
                loglike += expectTransitions[i][j] * log(pts[i][j]);
            }
        }
    }
    return loglike;
}

// write an objective function
double objective_function_asex(double * par)
{
    // (em is a global)
    const int n = em.hmm[0].n;
    assert(em.hmm[0].n == em.hmm[1].n);
    const int numHmmStates = em.numHmmStates;

    const double rho = par[0]*par[0];
    const double theta = par[1]*par[1];
    const double Td = par[2]*par[2];

    const int numParams = em.numFreeLambdas+3;

    double * lambdas = (double *)chmalloc(sizeof(double) * (n+1));

    int i, j;

    for(i = 0; i < em.lambdaCounts[0]; i++)
    {
        // set the first, fixed lambdas at 1
        lambdas[i] = 1.0;
    }
    int lambdaIdx = em.lambdaCounts[0];
    for(i = 0; i < em.numFreeLambdas; i++)
    {
        for(j = 0; j < em.lambdaCounts[i+1]; j++)
        {
            lambdas[lambdaIdx++] = par[i+3]*par[i+3];
        }
    }

    Hmm * scratchHmm = &(em.hmm[!em.hmmFlag]);

    Hmm_make_hmm(scratchHmm, lambdas, n, theta, rho, Td);

    double loglike = Em_get_expected_log_likelihood(&em, !em.hmmFlag);
    assert(loglike < 0);

    free(lambdas);

    return -loglike;
}

double objective_function_no_asex(double * par)
{
    // (em is a global)
    const int n = em.hmm[0].n;
    assert(em.hmm[0].n == em.hmm[1].n);
    const int numHmmStates = em.numHmmStates;

    const double rho = par[0]*par[0];
    const double theta = par[1]*par[1];

    const int numParams = em.numFreeLambdas+2;

    double * lambdas = (double *)chmalloc(sizeof(double) * (n+1));

    int i, j;

    for(i = 0; i < em.lambdaCounts[0]; i++)
    {
        // set the first, fixed lambdas at 1
        lambdas[i] = 1.0;
    }
    int lambdaIdx = em.lambdaCounts[0];
    for(i = 0; i < em.numFreeLambdas; i++)
    {
        for(j = 0; j < em.lambdaCounts[i+1]; j++)
        {
            lambdas[lambdaIdx++] = par[i+2]*par[i+2];
        }
    }

    Hmm * scratchHmm = &(em.hmm[!em.hmmFlag]);

    Hmm_make_hmm(scratchHmm, lambdas, n, theta, rho, 0.0);

    double loglike = Em_get_expected_log_likelihood(&em, !em.hmmFlag);
    assert(loglike < 0);

    free(lambdas);

    return -loglike;
}

double Em_get_loglikelihood(Em * em)
{
    double loglike = 0.0;
    const int numSeqs = em->numSeqs;
    const int numHmmStates = em->numHmmStates;
    int i, j, seqLen;
    double ** seqFor;
    double sum;
    for(i = 0; i < numSeqs; i++)
    {
        seqLen = em->dat->seqs[i].len;
        seqFor = em->forward[i];
        sum = 0.0;
        for(j = 0; j < numHmmStates; j++)
        {
            sum += seqFor[seqLen-1][j];
        }
        loglike += log(sum);

        for(j = 0; j < seqLen; j++)
        {
            loglike += log(em->normConst[i][j]);
        }
    }
    return loglike;
}

void Em_iterate(Em * em)
{
    if(em->asexEnabled)
    {
        Em_iterate_asex(em);
    }
    else
    {
        Em_iterate_no_asex(em);
    }
    return;
}

void Em_iterate_asex(Em * em)
{
    em->curIteration++;

    timestamp("starting forward-backward...");
    Em_get_forward(em);
    Em_get_backward(em);
    Em_get_expectations(em);
    timestamp("finished forward-backward...");

    assert(em->hmm[0].n == em->hmm[1].n);
    const int n = em->hmm[0].n;

    const int numOptimStarts = 20;

    double fmin;
    int i, j;

    int numParams = em->numFreeLambdas + 3;

    double * start = (double *)chmalloc(sizeof(double) * numParams);
    start[0] = sqrt(em->hmm[em->hmmFlag].rho);
    start[1] = sqrt(em->hmm[em->hmmFlag].theta);
    start[2] = sqrt(em->hmm[em->hmmFlag].Td);
    double * step = (double *)chmalloc(sizeof(double) * numParams);
    step[0] = sqrt(0.1);
    step[1] = sqrt(0.1);
    step[2] = sqrt(0.1);
    for(i = 3; i < numParams; i++)
    {
        start[i] = sqrt(em->freeLambdas[i-3]);
        step[i] = sqrt(0.5);
    }

    double ** randstarts = (double **)chmalloc(sizeof(double *) * numOptimStarts);
    randstarts[0] = (double *)chmalloc(sizeof(double) * numParams);
    for(j = 0; j < numParams; j++)
    {
        randstarts[0][j] = start[j];
    }
    double ** steps = (double **)chmalloc(sizeof(double *) * numOptimStarts);
    for(i = 0; i < numOptimStarts; i++)
    {
        steps[i] = (double *)chmalloc(sizeof(double) * numParams);
    }
    for(i = 0; i < numParams; i++)
    {
        steps[0][i] = step[i];
    }
    for(i = 1; i < numOptimStarts; i++)
    {
        for(j = 0; j < numParams; j++)
        {
            steps[i][j] = pow(2, runifab(-2, 2)) * steps[0][j];
        }
    }
    for(i = 1; i < numOptimStarts; i++)
    {
        randstarts[i] = (double *)chmalloc(sizeof(double) * numParams);
        for(j = 0; j < numParams; j++)
        {
            randstarts[i][j] = pow(10.0, runifab(-1.0, 1.0)) * start[j];
        }
    }

    int konvge = 1, maxNumEval = 10000;
    int iterationCount = 0, numRestarts, errorNum;
    double reqmin = 1e-16;

    double ** fargmins = (double **)chmalloc(sizeof(double *) * numOptimStarts);
    double * fmins = (double *)chmalloc(sizeof(double) * numOptimStarts);
    for(i = 0; i < numOptimStarts; i++)
    {
        timestamp("starting optimization rep...");
        fargmins[i] = (double *)chmalloc(sizeof(double) * numParams);
        nelmin(objective_function_asex, numParams, randstarts[i], fargmins[i],
                &(fmins[i]), reqmin, steps[i], konvge, maxNumEval, &iterationCount,
                &numRestarts, &errorNum);
        timestamp("optimization rep finished...");
    }

    for(i = 0; i < numOptimStarts; i++)
    {
        free(steps[i]);
    }
    free(steps);

    int minFminIdx = 0;
    double minfmin = fmins[0];
    for(i = 1; i < numOptimStarts; i++)
    {
        if(fmins[i] < minfmin)
        {
            minFminIdx = i;
            minfmin = fmins[i];
        }
    }

    double rho, theta, Td;

    rho = fargmins[minFminIdx][0]*fargmins[minFminIdx][0];
    theta = fargmins[minFminIdx][1]*fargmins[minFminIdx][1];
    Td = fargmins[minFminIdx][2]*fargmins[minFminIdx][2];

    for(i = 0; i < em->numFreeLambdas; i++)
    {
        em->freeLambdas[i] = fargmins[minFminIdx][i+3]*fargmins[minFminIdx][i+3];
    }

    double * lambdas = (double *)chmalloc(sizeof(double) * (n+1));

    for(i = 0; i < em->lambdaCounts[0]; i++)
    {
        // set the first, fixed lambdas at 1
        lambdas[i] = 1.0;
    }
    int lambdaIdx = em->lambdaCounts[0];
    for(i = 0; i < em->numFreeLambdas; i++)
    {
        for(j = 0; j < em->lambdaCounts[i+1]; j++)
        {
            lambdas[lambdaIdx++] = em->freeLambdas[i];
        }
    }
    assert(lambdaIdx == n+1);

    // set the HMM for the next iteration
    Hmm_make_hmm(&(em->hmm[!em->hmmFlag]), lambdas, n, theta, rho, Td);

    // flip the hmm flag
    em->hmmFlag = !em->hmmFlag;

    // free things
    free(lambdas);
    free(start);
    for(i = 0; i < numOptimStarts; i++)
    {
        free(fargmins[i]);
        free(randstarts[i]);
    }
    free(fargmins);
    free(randstarts);
    free(step);
    free(fmins);

    return;
}

void Em_iterate_no_asex(Em * em)
{
    em->curIteration++;

    //timestamp("making forward");
    Em_get_forward(em);
    //timestamp("making backward");
    Em_get_backward(em);
    //timestamp("making gamma and xi");
    Em_get_expectations(em);

    assert(em->hmm[0].n == em->hmm[1].n);
    const int n = em->hmm[0].n;

    const int numOptimStarts = 500;

    double fmin;
    int i, j;

    int numParams = em->numFreeLambdas + 2;

    double * start = (double *)chmalloc(sizeof(double) * numParams);
    start[0] = sqrt(em->hmm[em->hmmFlag].rho);
    start[1] = sqrt(em->hmm[em->hmmFlag].theta);
    double * step = (double *)chmalloc(sizeof(double) * numParams);
    step[0] = sqrt(0.1);
    step[1] = sqrt(0.1);
    for(i = 2; i < numParams; i++)
    {
        start[i] = sqrt(em->freeLambdas[i-2]);
        step[i] = sqrt(0.5);
    }

    double ** randstarts = (double **)chmalloc(sizeof(double *) * numOptimStarts);
    for(i = 0; i < numOptimStarts; i++)
    {
        randstarts[i] = (double *)chmalloc(sizeof(double) * numParams);
        for(j = 0; j < numParams; j++)
        {
            randstarts[i][j] = runifab(0.1, 2.0) * start[j];
        }
    }

    int konvge = 1, maxNumEval = 1000000;
    int iterationCount = 0, numRestarts, errorNum;
    double reqmin = 1e-10;

    double ** fargmins = (double **)chmalloc(sizeof(double *) * numOptimStarts);
    double * fmins = (double *)chmalloc(sizeof(double) * numOptimStarts);
    for(i = 0; i < numOptimStarts; i++)
    {
        fargmins[i] = (double *)chmalloc(sizeof(double) * numParams);
        nelmin(objective_function_no_asex, numParams, randstarts[i], fargmins[i],
                &(fmins[i]), reqmin, step, konvge, maxNumEval, &iterationCount,
                &numRestarts, &errorNum);
    }

    int minFminIdx = 0;
    double minfmin = fmins[0];
    for(i = 1; i < numOptimStarts; i++)
    {
        if(fmins[i] < minfmin)
        {
            minFminIdx = i;
            minfmin = fmins[i];
        }
    }

    em->hmm[!em->hmmFlag].rho = fargmins[minFminIdx][0]*fargmins[minFminIdx][0];
    em->hmm[!em->hmmFlag].theta = fargmins[minFminIdx][1]*fargmins[minFminIdx][1];

    for(i = 0; i < em->numFreeLambdas; i++)
    {
        em->freeLambdas[i] = fargmins[minFminIdx][i+2]*fargmins[minFminIdx][i+2];
    }

    for(i = 0; i < em->lambdaCounts[0]; i++)
    {
        // set the first, fixed lambdas at 1
        em->hmm[!em->hmmFlag].lambdas[i] = 1.0;
    }
    int lambdaIdx = em->lambdaCounts[0];
    for(i = 0; i < em->numFreeLambdas; i++)
    {
        for(j = 0; j < em->lambdaCounts[i+1]; j++)
        {
            em->hmm[!em->hmmFlag].lambdas[lambdaIdx++] = em->freeLambdas[i];
        }
    }

    em->hmmFlag = !em->hmmFlag;

    free(start);
    for(i = 0; i < numOptimStarts; i++)
    {
        free(fargmins[i]);
        free(randstarts[i]);
    }
    free(fargmins);
    free(randstarts);
    free(step);
    free(fmins);

    return;
}

void Em_print_forward(Em * em)
{
    const int numSeqs = em->numSeqs;
    const int numHmmStates = em->hmm[0].numStates;
    assert(em->hmm[0].numStates == em->hmm[1].numStates);
    int i, j, k, seqLen;
    double ** seqFor;
    Hmm * hmm = &(em->hmm[em->hmmFlag]);
    Data * dat = em->dat;

    for(i = 0; i < numSeqs; i++)
    {
        seqFor = em->forward[i];
        seqLen = dat->seqs[i].len;
        for(j = 0; j < seqLen; j++)
        {
            for(k = 0; k < numHmmStates; k++)
            {
                printf("%.15f\t", seqFor[j][k]);
            }
            printf("\n");
        }
    }
    return;
}

void Em_print_backward(Em * em)
{
    const int numSeqs = em->numSeqs;
    const int numHmmStates = em->hmm[0].numStates;
    assert(em->hmm[0].numStates == em->hmm[1].numStates);
    int i, j, k, seqLen;
    double ** seqBack;
    Hmm * hmm = &(em->hmm[em->hmmFlag]);
    Data * dat = em->dat;

    for(i = 0; i < numSeqs; i++)
    {
        seqBack = em->backward[i];
        seqLen = dat->seqs[i].len;
        for(j = 0; j < seqLen; j++)
        {
            for(k = 0; k < numHmmStates; k++)
            {
                printf("%.15f\t", seqBack[j][k]);
            }
            printf("\n");
        }
    }
    return;
}

void Em_print_norm_const(Em * em)
{
    const int numSeqs = em->numSeqs;
    const int numHmmStates = em->hmm[0].numStates;
    assert(em->hmm[0].numStates == em->hmm[1].numStates);
    int i, j, k, seqLen;
    double * seqNorm;
    Hmm * hmm = &(em->hmm[em->hmmFlag]);
    Data * dat = em->dat;

    for(i = 0; i < numSeqs; i++)
    {
        seqNorm = em->normConst[i];
        seqLen = dat->seqs[i].len;
        for(j = 0; j < seqLen; j++)
        {
            printf("%.15f\t", seqNorm[j]);
            printf("\n");
        }
    }
    return;
}

void Em_print_gamma(Em * em)
{
    const int numSeqs = em->numSeqs;
    const int numHmmStates = em->hmm[0].numStates;
    assert(em->hmm[0].numStates == em->hmm[1].numStates);
    int i, j, k, seqLen;
    double ** seqGamma;
    Hmm * hmm = &(em->hmm[em->hmmFlag]);
    Data * dat = em->dat;

    for(i = 0; i < numSeqs; i++)
    {
        seqGamma = em->gamma[i];
        seqLen = dat->seqs[i].len;
        for(j = 0; j < seqLen; j++)
        {
            for(k = 0; k < numHmmStates; k++)
            {
                printf("%.15f\t", seqGamma[j][k]);
            }
            printf("\n");
        }
    }
    return;
}

void Em_print_expect(Em * em)
{
    const int numSeqs = em->numSeqs;
    const int numHmmStates = em->hmm[0].numStates;
    assert(em->hmm[0].numStates == em->hmm[1].numStates);
    int i, j, k, seqLen;
    double ** seqGamma;
    Hmm * hmm = &(em->hmm[em->hmmFlag]);
    Data * dat = em->dat;

    for(i = 0; i < numHmmStates; i++)
    {
        for(j = 0; j < numHmmStates; j++)
        {
            printf("%.15f\t", em->expectTransitions[i][j]);
        }
        printf("\n");
    }
    return;
}

void Em_print_parameters(Em * em)
{
    Hmm * hmm = &(em->hmm[em->hmmFlag]);
    printf("PA\trho\t%g\n", hmm->rho);
    printf("PA\ttheta\t%g\n", hmm->theta);
    if(em->asexEnabled)
    {
        printf("PA\tTd\t%g\n", hmm->Td);
    }
    int i;
    for(i = 0; i < hmm->n+1; i++)
    {
        printf("PA\tlam%i\t%g\t%g\n", i, hmm->ts[i], hmm->lambdas[i]);
    }
    return;
}

void Em_print_iteration(Em * em)
{
    printf("IT\t%i\n", em->curIteration);
    printf("LL\t%g\n", Em_get_loglikelihood(em));
    Em_print_parameters(em);
    printf("CC\t---------------\n");
    return;
}

