#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "definitions.h"
#include "data.h"

/* data format:
 *   FASTA-like format, with each chromosome looking like
 *
 *   >comment for seq 1
 *   00000100002010001000000000001001002000000
 *   00010010020000000001002000000000000010
 *   >comment for seq 2
 *   00000100002010000000000000001001002000000
 *   00010010020000000001002000101200010
 *
 *   etc.
 *   
 */

const int DEFAULTINCREMENT = 1000;
const int NUMSEQINCREMENT = 100;

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

inline void Seq_add_buffer_size(Seq * seq, int increment)
{
    assert(seq);
    seq->data = (char *)chrealloc((void *)seq->data, sizeof(char) * (seq->maxLen + increment));
    seq->maxLen += increment;
}

inline void Seq_increment_length(Seq * seq)
{
    assert(seq);
    assert(seq->len < seq->maxLen);

    seq->len++;
    if(seq->len == seq->maxLen)
    {
        Seq_add_buffer_size(seq, DEFAULTINCREMENT);
    }
    return;
}

void Seq_init(Seq * seq, int initSize, int index)
{
    seq->data = (char *)chmalloc(sizeof(char) * initSize);
    seq->maxLen = initSize;
    seq->len = 0;
    seq->index = index;
    return;
}

/* Seq_read_seq
 *
 * takes a pointer to a FILE stream, reads a sequence of 0,1,[2,3], which
 * represent the data.
 *
 * sets fin to next > if it exists, or to EOF if reached
 *
 * seq    pointer to Seq, just initialized, maxLen must be 0
 * fin    pointer to stream, must be set to first character of data line
 * type   type of sequence, must be polarized or nonpolarized
 *
 */
void Seq_read_seq(Seq * seq, FILE * fin, SeqType type)
{
    assert(seq && fin);
    assert(seq->len == 0);

    assert(seq->len < seq->maxLen);

    const int maxValue = ((type == polarized) ? 3 : 1);

    int ch;
    while((ch = fgetc(fin)) != EOF)
    {
        if(ch == '\n')
        {
            continue;
        }
        if(ch == '>')
        {
            // rewind one so parent function sees '>'
            fseek(fin, -1, SEEK_CUR);
            return;
        }
        else
        {
            ch -= '0';
            assert(!(ch < 0 || ch > maxValue));
            if(ch < 0 || ch > maxValue)
            {
                PERROR("Invalid character in input");
            }
            seq->data[seq->len] = ch; // store the data
            // increment length of Seq
            Seq_increment_length(seq);
        }
    }
    return;
}

void Seq_print_seq(Seq * seq, int idx)
{
    assert(seq);
    int i;
    const int width = 60;
    printf(">%i\n", idx);
    assert(seq->len > 0);
    for(i = 0; i < seq->len; i++)
    {
        printf("%i", seq->data[i]);
        if(i > 0 && i % 60 == 0)
        {
            printf("\n");
        }
    }
    printf("\n");
    return;
}

void Seq_free(Seq * seq)
{
    assert(seq->data);
    free(seq->data);
    return;
}

void Data_init(Data * dat, SeqType seqtype)
{
    assert(dat);
    dat->numSeqs = 0;
    dat->seqs = (Seq *)chmalloc(sizeof(Seq) * NUMSEQINCREMENT);
    int i;
    dat->maxNumSeqs = NUMSEQINCREMENT;
    dat->seqtype = seqtype;
    return;
}

void Data_add_seq_buffer(Data * dat, int increment)
{
    assert(dat);
    assert(dat->numSeqs <= dat->maxNumSeqs);
    dat->seqs = chrealloc(dat->seqs, dat->maxNumSeqs + increment);
    dat->maxNumSeqs += increment;
    return;
}

// note that this function increments dat->numSeqs.
Seq * Data_get_seq(Data * dat)
{
    assert(dat);
    Seq * seq;
    if(dat->numSeqs == dat->maxNumSeqs)
    {
        Data_add_seq_buffer(dat, NUMSEQINCREMENT);
    }
    assert(dat->numSeqs < dat->maxNumSeqs);
    seq = &(dat->seqs[dat->numSeqs]);
    dat->numSeqs++;
    assert(seq);
    return seq;
}

void Data_read_seq(Data * dat, FILE * fin)
{
    assert(dat && fin);
    Seq * seq = Data_get_seq(dat);
    Seq_init(seq, DEFAULTINCREMENT, dat->numSeqs);
    Seq_read_seq(seq, fin, dat->seqtype);
}


/* Data_read_data
 *
 * reads sequence data from a file stream pointed to by fin, stores in dat
 *
 * dat   pointer to initialized Data object
 * fin   pointer to file stream containing data
 *
 */
void Data_read_data(Data * dat, FILE * fin)
{
    int ch, i, description;
    int chromLen, curChrom, status;
    long curPos;

    curChrom = -1;

    while((ch = getc(fin)) != EOF)
    {
        if(ch == '>')
        {
            // advance to end of line, loop will change to first character of
            // data in next iteration
            ch = finish_line(fin);
            Data_read_seq(dat, fin);
            continue;
        }
    }
    return;
}

void Data_print_seqs(Data * dat)
{
    int i;
    for(i = 0; i < dat->numSeqs; i++)
    {
        Seq_print_seq(&(dat->seqs[i]), i);
    }
    return;
}

