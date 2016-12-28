#ifndef GENOTYPELIKELIHOODS_H
#define GENOTYPELIKELIHOODS_H

/*
 * 
 * Data structure
 * ========================================
 *
 * GenotypeLikeData
 *
 *  - num_seqs
 *  - sequences (containing the genotype likelihoods)
 *  - current number of allocated slots
 *
 * SeqGenotypeLikelihood
 *  - length
 *  - matrix of log-likelihoods of the emission genotypes at that site
 *      . for non-polarized data, just have total log-like of homozygous and
 *        heterozygous
 *      . for polarized data, have log-likelihoods of homozygous, or 1 or 2
 *        ancestral copies
 *  - name
 *   
 * */

struct genotypelikedata_;

typedef struct seqgenotypelikelihood_
{
    int len;
    int max_len;
    int num_emission_states;
    double *likes;  // likes is len x num_emission_states in length
    struct genotypelikedata_ *parent;
    char * name; // sequence name
} SeqGenotypeLike;

typedef struct genotypelikedata_
{
    int num_seqs;
    int num_allocated_slots;
    SeqGenotypeLike *seq_likes;
    int num_emission_states;
    int polarized;
} GenotypeLikeData;


void GenotypeLikeData_init(GenotypeLikeData *dat, int polarized);
void GenotypeLikeData_free(GenotypeLikeData *dat);
void SeqGenotypeLike_free(SeqGenotypeLike *seq);
void SeqGenotypeLike_add_likelihoods(SeqGenotypeLike *seq,
    double *m_linelikes);
SeqGenotypeLike * GenotypeLikeData_get_next_seq(GenotypeLikeData *dat,
        char *name);
void GenotypeLikeData_read_data(FILE *fin, GenotypeLikeData *dat,
        int polarized);
void SeqGenotypeLike_init(SeqGenotypeLike *seq,
        char *name, GenotypeLikeData *parent);
void SeqGenotypeLike_free(SeqGenotypeLike *seq);
void SeqGenotypeLike_add_likelihoods(SeqGenotypeLike *seq,
    double *m_linelikes);
SeqGenotypeLike * GenotypeLikeData_get_next_seq(GenotypeLikeData *dat,
        char *name);
void GenotypeLikeData_read_data(FILE *fin, GenotypeLikeData *dat,
        int polarized);
void GenotypeLikeData_print(GenotypeLikeData *dat, FILE *fout);


#endif
