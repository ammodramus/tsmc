#include <stdio.h>
#include <stdlib.h>
#include "definitions.h"

/* data format:
 *   FASTA-like format, with each chromosome looking like
 *
 *   >chromosome description
 *   00000100002010001000000000001001002000000
 *   00010010020000000001002000000000000000010
 *
 *   etc.
 *   
 */

int finish_line(FILE * fin)
{
    int ch;
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

void Data_allocate_chromosome(Data * dat, int chromLen, int chromIdx)
{
    return;
}

void Data_read_data(FILE * fin, Data * dat)
{
    int ch, i, description;
    int chromLen, curChrom, status;
    long curPos;

    curChrom = -1;

    while((ch = getc(fin)) != EOF)
    {
        if(ch == '>')
        {
            curChrom++;

            ch = finish_line(fin);

            // now read chromosome length
            curPos = ftell(fin);
            chromLen = get_chrom_len(fin);
            fseek(fin, curPos, SEEK_SET);

            // allocate space for current chromosome
            Data_allocate_chromosome(dat, chromLen, curChrom);
            
            //rewind and store data

            continue;
        }
        else if(ch == '\n')
        {
            continue;
        }
        else
        {


    return;
}

