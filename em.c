#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "definitions.h"
#include "hmm.h"
#include "data.h"
#include "em.h"

void Em_init(Em * em, Data * dat, Hmm * hmm)
{
    assert(em && dat && hmm);
    assert(dat->numSeqs > 0);
    em->dat = dat;
    em->hmm = hmm;
    em->numSeqs = dat->numSeqs;
    em->seqtype = dat->seqtype;
    em->forward = (double ***)chmalloc(sizeof(double **) * em->numSeqs);
    em->backward = (double ***)chmalloc(sizeof(double **) * em->numSeqs);
    em->normConst = (double **)chmalloc(sizeof(double *) * em->numSeqs);
    em->gamma = (double ***)chmalloc(sizeof(double **) * em->numSeqs);
    int i, j;
    for(i = 0; i < em->numSeqs; i++)
    {
        em->forward[i] = (double **)chmalloc(sizeof(double *) * dat->seqs[i].len);
        em->backward[i] = (double **)chmalloc(sizeof(double *) * dat->seqs[i].len);
        em->gamma[i] = (double **)chmalloc(sizeof(double *) * dat->seqs[i].len);
        em->normConst[i] = (double *)chmalloc(sizeof(double) * dat->seqs[i].len);
        for(j = 0; j < dat->seqs[i].len; j++)
        {
            em->forward[i][j] = (double *)chmalloc(sizeof(double) * em->hmm->numStates);
            em->backward[i][j] = (double *)chmalloc(sizeof(double) * em->hmm->numStates);
            em->gamma[i][j] = (double *)chmalloc(sizeof(double) * em->hmm->numStates);
        }
    }
    em->expect = (double **)chmalloc(sizeof(double *) * hmm->numStates);
    for(i = 0; i < hmm->numStates; i++)
    {
        em->expect[i] = (double *)chmalloc(sizeof(double) * hmm->numStates);
    }
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
    }
    free(em->forward);
    free(em->backward);
    free(em->gamma);
    return;
}

void Em_get_forward(Em * em)
{
    double * const pis = em->hmm->pis;
    double ** const pts = em->hmm->pts;
    double *** const forward = em->forward;
    const int numEmissionStates = (em->seqtype == polarized) ? 4 : 2;
    const int numHmmStates = em->hmm->numStates;
    const int numSeqs = em->numSeqs;
    fourd * const emissions = em->hmm->emissions;

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
        sum1 = 0.0;
        for(k = 0; k < numHmmStates; k++)
        {
            seqFor[0][k] = pis[k]*emissions[k][seqData[0]];
            sum1 += seqFor[0][k];
        }
        em->normConst[i][0] = sum1;
        for(k = 0; k < numHmmStates; k++)
        {
            seqFor[0][k] /= sum1;
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
                    sum1 += seqFor[j-1][l]*pts[l][k];
                }
                seqFor[j][k] = emissions[k][seqData[j]] * sum1;
                sum2 += seqFor[j][k];
            }
            em->normConst[i][j] = sum2;
            for(k = 0; k < numHmmStates; k++)
            {
                seqFor[j][k] /= sum2;
            }
        }
    }
    return;
}

void Em_get_backward(Em * em)
{
    double * const pis = em->hmm->pis;
    double ** const pts = em->hmm->pts;
    double *** const backward = em->backward;
    const int numEmissionStates = (em->seqtype == polarized) ? 4 : 2;
    const int numHmmStates = em->hmm->numStates;
    const int numSeqs = em->numSeqs;
    fourd * const emissions = em->hmm->emissions;

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
                    sum += seqBack[j+1][l]*pts[k][l]*emissions[l][seqData[j+1]];
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
    double *** const forward = em->forward;
    double *** const backward = em->backward;
    double *** const gamma = em->gamma;
    double ** const expect = em->expect;
    double ** const pts = em->hmm->pts;
    fourd * const emissions = em->hmm->emissions;
    const int numHmmStates = em->hmm->numStates;
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
            expect[j][k] = 0.0;
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
                sum += seqFor[j][l]*seqBack[j][l];
            }
            for(l = 0; l < numHmmStates; l++)
            {
                em->gamma[i][j][l] /= sum;
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
                    thisExpect = em->gamma[i][j][k]*pts[k][l]*seqBack[j+1][l] *
                        emissions[l][seqData[j+1]] / (em->normConst[i][j+1] * seqBack[j][k]);
                    expect[k][l] += thisExpect;
                    assert(0 <= thisExpect && thisExpect <= 1);
                }
            }
        }
    }
    return;
}
