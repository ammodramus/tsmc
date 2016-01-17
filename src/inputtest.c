#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "hmm.h"
#include "definitions.h"
#include "data.h"
#include "em.h"

int main(int argc, char ** argv)
{
    Data dat;
    Data_init(&dat, polarized);
    FILE * fin;
    if(argc == 1)
    {
        fin = chfopen("testseqs", "r");
    }
    else if(argc == 2)
    {
        fin = chfopen(argv[1], "r");
    }
    Data_read_data(&dat, fin);
    printf("numSeqs = %i\n", dat.numSeqs);
    Data_print_seqs(&dat);
    Data_free(&dat);
    fclose(fin);
    return 0;
}

