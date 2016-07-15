#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#ifdef DEBUGMEMORY
#include <malloc.h>
#endif
#include "definitions.h"


#ifndef DEBUGMEMORY
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

inline void chfree(void * ptr)
{
    if(ptr == NULL)
    {
        fprintf(stderr, "Double free attempted.\n");
        exit(1);
    }
    free(ptr);
    ptr = NULL;
}
#endif

#ifdef DEBUGMEMORY
inline void * chmalloc_memdebug(size_t size, char * filename, int lineNumber)
{
    void * ptr = malloc(size);
    fprintf(stdout, "all %i %s %i\n", (int)size, filename, lineNumber);
    if(!ptr)
    {
        fprintf(stderr, "Out of memory.\n");
        exit(1);
    }
    return ptr;
}

inline void * chrealloc_memdebug(void * oldptr, size_t size, char * filename, int lineNumber)
{
    void * ptr = realloc(oldptr, size);
    fprintf(stdout, "reall %i %s %i\n", (int)size, filename, lineNumber);
    if(!ptr)
    {
        fprintf(stderr, "Out of memory.\n");
        exit(1);
    }
    return ptr;
}

inline void * chcalloc_memdebug(size_t nmemb, size_t size, char * filename, int lineNumber)
{
    void * ptr = calloc(nmemb, size);
    fprintf(stdout, "call %i %s %i\n", (int)size, filename, lineNumber);
    if(!ptr)
    {
        fprintf(stderr, "Out of memory.\n");
        exit(1);
    }
    return ptr;
}

inline void * chfree_memdebug(void *ptr, char * filename, int lineNumber)
{
    if(ptr == NULL)
    {
        fprintf(stderr, "Double free attempted.\n");
        exit(1);
    }
    fprintf(stdout, "free %i %s %i\n", (int)malloc_usable_size(ptr), filename, lineNumber);
    free(ptr);
}

#endif


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

#ifndef NDEBUG
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
#else
inline void timestamp(const char * msg)
{
    return;
}
#endif

void get_ts_msmc(double * ts, const int n)
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

void get_ts_psmc(double * ts, const double tmax, const double n)
{
    int i;
    ts[0] = 0.0;
    for(i = 1; i <= n; i++)
    {
        ts[i] = 0.1*(exp((double)i/n*log(1.0+10.0*tmax))-1.0);
    }
    return;
}

