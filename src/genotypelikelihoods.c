#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "definitions.h"
#include "genotypelikelihoods.h"
#include "data.h"

/* data format:
 * contig   pos likegen0    likeg1  ...     likeg19
 */

static const int DEFAULTINCREMENT = 1000;
static const int NUMSEQINCREMENT = 100;
static const int INITIALSIZE = 1000;
static const int DEFAULTINCREMENTFACTOR = 2;

static const int MAX_NAME = 1024;

int finish_line(FILE * fin)
{
    int ch;
    ch = fgetc(fin);
    while(ch != '\n')
    {
        ch = fgetc(fin);
        if(ch == EOF)
        {
            return(EOF);
        }
    }
    return ch;
}

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

void merge_line_likes(double *linelikes, double *m_linelikes, int polarized)
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
