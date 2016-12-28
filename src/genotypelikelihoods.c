#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "definitions.h"
#include "genotypelikelihoods.h"
#include "data.h"
#include "util.h"
#include "hmm.h"

/* data format:
 * contig   pos likeg0    likeg1  ...     likeg19
 */

static const int DEFAULTINCREMENT = 1000;
static const int NUMSEQINCREMENT = 100;
static const int INITIALSIZE = 1000;
static const int DEFAULTINCREMENTFACTOR = 2;

static const int MAX_NAME = 1024;

void GenotypeLikeData_init(GenotypeLikeData *dat, int polarized)
{
    dat->polarized = polarized;
    if (polarized)
    {
        dat->num_emission_states = 3;
    }
    else
    {
        dat->num_emission_states = 2;
    }
    dat->seq_likes = chmalloc(sizeof(SeqGenotypeLike) * INITIALSIZE);
    dat->num_allocated_slots = INITIALSIZE;
    dat->num_seqs = 0;
    return;
}

void GenotypeLikeData_free(GenotypeLikeData *dat)
{
    int i;
    for (i = 0; i < dat->num_seqs; i++)
    {
        SeqGenotypeLike_free(&(dat->seq_likes[i]));
    }
    free(dat->seq_likes);
    return;
}

static int HETEROZYGOUS[] = {
    // AAA AAC AAG AAT ACC ACG ACT AGG AGT ATT CCC CCG CCT CGG CGT CTT GGG GGT
       0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  1,  1,  1,  1,  1,  0,  1,
    // GTT TTT
       1,  0};

static void merge_line_likes(double *linelikes, double *m_linelikes, int polarized)
{
    int i;
    if (polarized)
    {
        // if polarized, need to know the ancestral base
        PERROR("polarized not implemented");
    }
    else
    {
        // if not polarized, can just use heterozygous vs. homozygous
        m_linelikes[0] = 0.0;
        m_linelikes[1] = 0.0;
        for (i = 0; i < 20; i++) 
        {
            m_linelikes[HETEROZYGOUS[i]] += exp(linelikes[i]);
        }
        // assume that P(G | het) = 1/16 for heterozygous genotypes
        m_linelikes[0] = log(m_linelikes[0] / 16.0);
        // assume that P(G | het) = 1/4 for homozygous genotypes
        m_linelikes[1] = log(m_linelikes[1] / 4.0);
        // but what about GC content, e.g.?
    }
}

void SeqGenotypeLike_init(SeqGenotypeLike *seq,
        char *name, GenotypeLikeData *parent)
{
    seq->len = 0;
    seq->max_len = INITIALSIZE;
    seq->likes = chcalloc(seq->max_len * parent->num_emission_states,
            sizeof(double));
    seq->name = chcalloc(MAX_NAME, sizeof(char));
    seq->name = strncpy(seq->name, name, MAX_NAME);
    seq->parent = parent;
    seq->num_emission_states = parent->num_emission_states;
    return;
}

void SeqGenotypeLike_free(SeqGenotypeLike *seq)
{
    free(seq->likes);
    free(seq->name);
    return;
}

void SeqGenotypeLike_add_likelihoods(SeqGenotypeLike *seq,
    double *m_linelikes)
{
    int i;
    for (i = 0; i < seq->num_emission_states; i++)
    {
        seq->likes[seq->num_emission_states*seq->len+i] = m_linelikes[i];
    }
    seq->len++;
    if (seq->len >= seq->max_len)
    {
        int new_len = (int)(seq->max_len*DEFAULTINCREMENTFACTOR);
        seq->likes = chrealloc(seq->likes,
                new_len * seq->num_emission_states * sizeof(double));
        seq->max_len = new_len;
    }
    return;
}

void SeqGenotypeLike_print(SeqGenotypeLike *seq, FILE *fout)
{
    int i, j;
    for (i = 0; i < seq->len; i++)
    {
        fprintf(fout, "%s\t", seq->name);
        fprintf(fout, "%i\t", i);
        for (j = 0; j < seq->num_emission_states; j++)
        {
            fprintf(fout, "%f\t", seq->likes[i*seq->num_emission_states+j]);
        }
        fprintf(fout, "\n");
    }
    return;
}

SeqGenotypeLike *GenotypeLikeData_get_next_seq(GenotypeLikeData *dat,
        char *name)
{
    SeqGenotypeLike *next_seq;
    if (dat->num_allocated_slots <= dat->num_seqs+1)
    {
        int new_num_allocations = (int)(dat->num_allocated_slots *
                DEFAULTINCREMENTFACTOR);
        dat->seq_likes = chrealloc(dat->seq_likes,
                sizeof(SeqGenotypeLike));
        dat->num_allocated_slots = new_num_allocations;
    }
    next_seq = &(dat->seq_likes[dat->num_seqs]);
    SeqGenotypeLike_init(next_seq, name, dat);
    dat->num_seqs++;
    return next_seq;
}
    

void GenotypeLikeData_read_data(FILE *fin, GenotypeLikeData *dat,
        int polarized)
{
    char line[MAX_NAME], namebuf[MAX_NAME], prev_name_buf[MAX_NAME];
    char *name, *prev_name, *tmp_line, *tok;
    int num_likes, pos, i;
    double *linelikes, *m_linelikes;
    SeqGenotypeLike *seq;

    num_likes = 20;
    linelikes = chcalloc(num_likes, sizeof(double));
    if (polarized)
    {
        // merged line log-likelihoods, 0, 1, 2
        m_linelikes = chcalloc(3, sizeof(double));
    }
    else
    {
        // two states if not polarized (hom or het)
        m_linelikes = chcalloc(2, sizeof(double));
    }

    dat->polarized = polarized;

    // header
    fgets(line, MAX_NAME, fin);

    prev_name = NULL;
    prev_name = strncpy(prev_name_buf, "asdfasdf", MAX_NAME);
    while (fgets(line, MAX_NAME, fin))
    {
        tmp_line = strdup(line);
        tok = strtok(tmp_line, "\t");
        name = strncpy(namebuf, tok, MAX_NAME);
        if (strncmp(name, prev_name, MAX_NAME) != 0)
        {
            seq = GenotypeLikeData_get_next_seq(dat, name);
        }
        tok = strtok(NULL, "\t");
        pos = (int)strtol(tok, NULL, 10);
        for (i = 0; i < num_likes; i++)
        {
            tok = strtok(NULL, "\t");
            linelikes[i] = strtod(tok, NULL);
        }
        merge_line_likes(linelikes, m_linelikes, polarized);
        SeqGenotypeLike_add_likelihoods(seq, m_linelikes);
        prev_name = strncpy(prev_name_buf, name, MAX_NAME);
        free(tmp_line);
    }
out:
    free(linelikes);
    free(m_linelikes);
    return;
}

void GenotypeLikeData_print(GenotypeLikeData *dat, FILE *fout)
{
    int i;
    for (i = 0; i < dat->num_seqs; i++)
    {
        SeqGenotypeLike_print(&(dat->seq_likes[i]), fout);
    }
    return;
}

/*
 * SeqGenLikeEmissions
 *
 */

void SeqGenLikeEmissions_init(SeqGenLikeEmissions *se, int num_bins,
        int num_hidden_states, int bin_width);
{
    int i;


    se->num_hidden_states = num_hidden_states;
    se->num_bins = num_bins;
    se->num_emission_states = num_emission_states;
    se->loglikes = chcalloc(num_bins, sizeof(double *));
    se->bin_width = bin_width;
    //
    // se->loglikes[i][j] will be P(S_i | T = j), where S_i is the reads in the
    // ith window of the binned sequence and j is the hidden state (i.e.,
    // ijidx).
    // 
    for (i = 0; i < num_bins; i++)
    {
        se->loglikes[i] = chcalloc(num_hidden_states, sizeof(double));
    }
    return;
}

void SeqGenLikeEmissions_free(SeqGenLikeEmissions *se)
{
    int i;
    for (i = 0; i < se->num_bins; i++)
    {
        chfree(se->loglikes[i])
    }
    chfree(se->loglikes);
    return;
}

void SeqGenLikeEmissions_calculate(SeqGenLikeEmissions *se,
        SeqGenotypeLike *sl, double **obs_emissions, int num_hidden_states)
{
    int i, j, k, l, sl_pos, polarized, done;
    double bin_prob, site_prob;
    const int bin_width = se->bin_width;
    const int num_emission_states = se->num_emission_states;
    polarized = (se->num_emission_states == 3) ? 1 : 0;
    assert(se->num_emission_states == 3 || se->num_emission_states == 2);

    if(!polarized)
    {
        sl_pos = 0;
        // calculate an emission probability for each bin
        for(i = 0; i < se->num_bins; i++)
        {
            // ... for each hidden state
            for(k = 0; k < num_hidden_states; k++)
            {
                bin_loglike = 0.0;
                done = 0;
                // for each site in the bin
                for(j = 0; j < bin_width && !done; j++)
                {
                    site_prob = 0.0;
                    // for each genotype...
                    for(l = 0; l < num_emission_states; l++)
                    {
                        site_prob +=
                            exp(sl->likes[sl_pos*num_emission_states + l]) *
                            obs_emissions[k][l];
                    }
                    bin_loglike += log(site_prob);
                    sl_pos++;
                    if(sl_pos >= sl->len)
                    {
                        done = 1;
                        break;
                    }
                }
                se->loglikes[i][k] = bin_loglike;
            }
        }
    }
    else
    {
        PERROR("polarized not yet implemented.");
    }
}



/*
 ******************************************************************
 *
 * GenLikeEmissions -- object for likelihoods that have been condensed
 *
 * Gives P(S^{(\nu)} | T) for each binned window \nu and each hidden state T.
 *
 *****************************************************************
 */


void GenLikeEmissions_calculate(GenLikeEmissions *emissions,
        GenotypeLikeData *dat, int bin_width)
{
    emissions->num_seqs = dat->num_seqs;

}
void GenLikeEmissions_calculate(GenLikeEmissions *emissions,
        GenotypeLikeData *dat, int bin_width, Hmm *hmm)
{
    int i, j, ijidx, num_bins;
    double Es3, Es2, tree_size;

    const int n = hmm->n;
    const int num_states = hmm->numStates;
    const double Td = hmm->Td;
    const double theta = hmm->theta;

    double *const Eijs2s = hmm->Eijs2s;
    double *const Eijs3s = hmm->Eijs3s;

    emissions->num_seqs = dat->num_seqs;
    emissions->bin_width = bin_width;
    emissions->num_emission_states = dat->num_emission_states;
    emissions->polarized = dat->polarized;
    emissions->seq_likes = chcalloc(dat->num_seqs,
            sizeof(SeqGenLikeEmissions));

    double **obs_emissions;

    obs_emissions = chcalloc(num_states, sizeof(double *));
    for(i = 0; i < num_states; i++)
    {
        obs_emissions[i] = chcalloc(6, sizeof(double));
    }


    for(i = 0; i <= n; i++)
    {
        for(j = i; j <= n; j++)
        {
            ijidx = get_index(i, j, n);
            Es2 = Eijs2s[ijidx];
            Es3 = Eijs3s [ijidx];
            tree_size = 2.0*Es2 + Es3;

            /* homozygous */
            obs_emissions[ijidx][0] = exp(-(tree_size + 3.0*Td)*theta/2.0);
            /* 1-mut */
            obs_emissions[ijidx][1] = exp(-(Es2-Es3)*theta/2.0) *
                (1 - exp(-(2*Es3 + Es2 + 3*Td)*theta/2.0));
            /* 2-mut */
            obs_emissions[ijidx][2] = (1.0 - exp(-(Es2-Es3)*theta/2.0)) * 
                exp(-(2*Es3 + Es2 + 3*Td)*theta/2.0);
            /* 1-mut and 2-mut */
            obs_emissions[ijidx][3] = (1.0 - exp(-(Es2-Es3)*theta/2.0)) * 
                (1 - exp(-(2*Es3 + Es2 + 3*Td)*theta/2.0));
            /* heterozygous, no polarization */
            obs_emissions[ijidx][4] = 1.0 - obs_emissions[ijidx][0];
            /* missing */
            obs_emissions[ijidx][5] = 1.0;
        }
    }

    if(!dat->polarized)
    {
        for(i = 0; i < dat->num_seqs; i++)
        {
            num_bins = dat->seq_likes[i].len / bin_width;
            SeqGenLikeEmissions_init(&(emissions->seq_likes[i]), num_bins,
                    num_states, bin_width);
            SeqGenLikeEmissions_calculate(&(emissions->seq_likes[i]),
                    &(dat->seq_likes[i]), obs_emissions, num_states);
        }
    }
    else
    {
        PERROR("polarized GenLikeEmissions not yet implemented.");
    }
    

out:
    for(i = 0; i < num_states; i++)
    {
        chfree(obs_emissions[i]);
    }
    chfree(obs_emissions);
    return;

}

void GenLikeEmissions_free(GenLikeEmissions *emissions)
{
    int i;
    for(i = 0; i < emissions->num_seqs; i++)
    {
        SeqGenLikeEmissions_free(&(emissions->seq_likes[i]));
    }
    chfree(emissions->seq_likes);
    return;
}
