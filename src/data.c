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
const int INITIALSIZE = 1000;
const double DEFAULTINCREMENTFACTOR = 1.4;

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

extern inline void Seq_add_buffer_size(Seq * seq, int increment)
{
    assert(seq);
    seq->data = (char *)chrealloc((void *)seq->data, sizeof(char) * (size_t)(seq->maxLen + increment));
    seq->maxLen += increment;
}

extern inline void Seq_multiply_buffer_size(Seq * seq, double factor)
{
    assert(seq);
    assert(factor > 1.0); // should be checked and not asserted...
    size_t nextLen = (size_t)((double)(seq->maxLen) * factor);
    seq->data = (char *)chrealloc((void *)seq->data, sizeof(char) * nextLen);
    seq->maxLen = nextLen;
}
extern inline void Seq_increment_length(Seq * seq)
{
    assert(seq);
    assert(seq->len < seq->maxLen);

    seq->len++;
    if(seq->len == seq->maxLen)
    {
        //Seq_add_buffer_size(seq, DEFAULTINCREMENT);
        Seq_multiply_buffer_size(seq, DEFAULTINCREMENTFACTOR);
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
int Seq_read_seq(Seq * seq, FILE * fin, SeqType type)
{
    assert(seq && fin);
    assert(seq->len == 0);

    assert(seq->len < seq->maxLen);

    /* allowed inputs are 0-4, and N */
    int foundHet;
    const int maxValue = 4;
    char ch;

    foundHet = 0;
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
            return foundHet;
        }
        if(ch == 'N')
        {
            /* N encodes for missing data in input, 5 in data */
            seq->data[seq->len] = 5;
            Seq_increment_length(seq);
        }
        else
        {
            ch -= '0';
            if(ch < 0 || ch > maxValue)
            {
                PERROR("Invalid character in input. Allowed characters are [01234N]");
            }
            seq->data[seq->len] = (char)ch; // store the data
            // increment length of Seq
            Seq_increment_length(seq);
            if(ch > 0)
            {
                foundHet = 1;
            }
        }
    }
    return foundHet;
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
        if(i > 0 && (i+1) % width == 0 && i < seq->len-1)
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
    chfree(seq->data);
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

void Data_free(Data * dat)
{
    int i;
    for(i = 0; i < dat->numSeqs; i++)
    {
        Seq_free(&(dat->seqs[i]));
    }
    chfree(dat->seqs);
    return;
}

void Data_add_seq_buffer(Data * dat, int increment)
{
    assert(dat);
    assert(dat->numSeqs <= dat->maxNumSeqs);
    dat->seqs = chrealloc(dat->seqs, sizeof(Seq) * (dat->maxNumSeqs + increment));
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

int Data_read_seq(Data * dat, FILE * fin)
{
    int foundHet;
    assert(dat && fin);
    Seq * seq = Data_get_seq(dat);
    Seq_init(seq, DEFAULTINCREMENT, dat->numSeqs);
    foundHet = Seq_read_seq(seq, fin, dat->seqtype);
    return foundHet;
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
    int ch, i, description, chromLen, curChrom, status, foundHet, seqFoundHet;
    long curPos;

    curChrom = -1;

    foundHet = 0;

    while((ch = getc(fin)) != EOF)
    {
        if(ch == '>')
        {
            // advance to end of line, loop will change to first character of
            // data in next iteration
            ch = finish_line(fin);
            seqFoundHet = Data_read_seq(dat, fin);
            if(seqFoundHet)
            {
                foundHet = 1;
            }
            continue;
        }
    }
    if(!foundHet)
    {
        PERROR("No heterozygous sites found in input data. At least one heterozygous site is required.");
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

