#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
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

inline void perror(const char * msg)
{
    fprintf(stderr, "%s\n", msg);
    exit(1);
}

void timestamp(const char * msg)
{
    char timestmp[100];
    struct tm * timeInfo;
    time_t rawTime;
    time(&rawTime);
    timeInfo = localtime(&rawTime);
    strftime(timestmp, 100, "%c", timeInfo);
    fprintf(stderr, "%s -- %s\n", timestmp, msg);
    return;
}

void get_ts(double * ts, const int n)
{
    double F;
    int i;
    ts[0] = 0.0;
    for(i = 1; i < n+1; i++)
    {
        F = (double)(i) * 1.0/((double)n+1.0);
        ts[i] = -log(1-F);
    }
    return;
}

