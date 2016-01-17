#ifndef DATA_H
#define DATA_H

#include "definitions.h"

/* Notes on data organization
 * ==========================
 *
 * 
 * Data struct needs:
 * -------------------
 * - number of seqs (int)
 * - type of data (polarized or not) (int or typedef'ed enum)
 * - data (Seq *)
 * - current number of allocated slots (int)
 *
 * Data methods/functions:
 * -----------------------
 *   - read data
 *   - allocate space for storing chroms
 *   - add a Seq
 *   - init, allocate a new Seq
 * 
 * Seq struct needs:
 * -----------------
 * - length (int)
 * - maxlength (int)
 * - index (int)
 * - data (char *)
 *
 * Seq methods/functions:
 * ----------------------
 * - read seq from FILE *
 * - increment Seq length
 * - (re)allocate space for seq
 */

// global increment size
typedef enum
{
    polarized = 0,
    nonpolarized
} SeqType;

typedef struct
{
    int len;
    int maxLen;
    int index;
    char * data;
} Seq;


typedef struct
{
    int numSeqs;
    int maxNumSeqs;
    SeqType seqtype;
    Seq * seqs;
} Data;

/*
 * Data struct needs:
 * -------------------
 * - number of seqs (int)
 * - type of data (polarized or not) (int or typedef'ed enum)
 * - data (Seq *)
 * - current number of allocated slots (int)
 *
 * Data methods/functions:
 * -----------------------
 *   - init
 *   - read data
 *   - allocate space for storing chroms
 *   - add a Seq
 *   - init, allocate a new Seq
 */

inline void Seq_add_buffer_size(Seq * seq, int increment);
inline void Seq_increment_length(Seq * seq);
void Seq_init(Seq * seq, int initSize, int index);
void Seq_read_seq(Seq * seq, FILE * fin, SeqType type);
void Seq_print_seq(Seq * seq, int idx);
void Seq_free(Seq * seq);
void Data_init(Data * dat, SeqType seqtype);
void Data_add_seq_buffer(Data * dat, int increment);
void Data_add_seq_buffer(Data * dat, int increment);
Seq * Data_get_seq(Data * dat);
void Data_read_seq(Data * dat, FILE * fin);
void Data_read_data(Data * dat, FILE * fin);
void Data_print_seqs(Data * dat);

#endif
