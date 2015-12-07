#include <stdio.h>
#include <stdlib.h>
#include "definitions.h"

void * chmalloc(size_t size)
{
    void * ptr = malloc(size);
    if(!ptr)
    {
        fprintf(stderr, "Out of memory.\n");
        exit(1);
    }
    return ptr;
}
