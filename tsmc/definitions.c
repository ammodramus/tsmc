#include <stdio.h>
#include <stdlib.h>
#include "definitions.h"

inline void * chmalloc(size_t size)
{
    void * ptr = malloc(size);
    if(!ptr)
    {
        fprintf(stderr, "Out of memory.\n");
        exit(1);
    }
    return ptr;
}

inline void * chrealloc(void * oldptr, size_t size)
{
    void * ptr = realloc(oldptr, size);
    if(!ptr)
    {
        fprintf(stderr, "Out of memory.\n");
        exit(1);
    }
    return ptr;
}

inline void * chcalloc(size_t nmemb, size_t size)
{
    void * ptr = calloc(nmemb, size);
    if(!ptr)
    {
        fprintf(stderr, "Out of memory.\n");
        exit(1);
    }
    return ptr;
}

inline FILE * chfopen(const char * path, const char * mode)
{
    FILE * fi = fopen(path, mode);
    if(!fi)
    {
        fprintf(stderr, "Could not open %s\n", path);
        exit(2);
    }
    return fi;
}
