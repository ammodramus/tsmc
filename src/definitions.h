#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <stdlib.h>

#define DEBUGREPORTI(x) do {fprintf(stderr, #x " = %i\n", x);} while(0)
#define DEBUGREPORTF(x) do {fprintf(stderr, #x " = %.15f\n", x);} while(0)
#define PERROR(msg,...) do {fprintf(stderr,"\n\nProgram error:\n%s\n\n",msg); exit(-1);} while(0)

typedef double fourd[4];

void * chmalloc(size_t size);
inline void * chrealloc(void * oldptr, size_t size);
inline void * chcalloc(size_t nmemb, size_t size);
inline FILE * chfopen(const char * path, const char * mode);
void timestamp(const char * msg);
void get_ts_msmc(double * ts, const int n);
void get_ts_psmc(double * ts, const double tmax, const double n);

#endif