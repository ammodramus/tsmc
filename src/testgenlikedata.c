#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "definitions.h"
#include "genotypelikelihoods.h"

//void GenotypeLikeData_read_data(FILE *fin, GenotypeLikeData *dat,
//        int polarized)

int main()
{
    int polarized;
    FILE *fin;
    polarized = 0;
    GenotypeLikeData dat;
    GenotypeLikeData_init(&dat, polarized);

    polarized = 0;
    fin = fopen("test_likelihoods.txt", "r");
    GenotypeLikeData_read_data(fin, &dat, 0);
    GenotypeLikeData_print(&dat, stdout);
    GenotypeLikeData_free(&dat);
    fclose(fin);

    return 0;
}
