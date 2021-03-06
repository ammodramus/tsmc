#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <stdlib.h>

#ifndef NDEBUG
#define DEBUGREPORTI(x) do {fprintf(stderr, #x " = %i\n", x);} while(0)
#define DEBUGREPORTF(x) do {fprintf(stderr, #x " = %.15f\n", x);} while(0)
#else
#define DEBUGREPORTI(x)
#define DEBUGREPORTF(x)
#endif
#define PERROR(msg,...) do {fprintf(stderr,"\n\nProgram error:\n%s\n\n",msg); exit(-1);} while(0)

typedef double sixd[6];

#ifndef DEBUGMEMORY
void * chmalloc(size_t size);
inline void * chrealloc(void * oldptr, size_t size);
inline void * chcalloc(size_t nmemb, size_t size);
inline void chfree(void *ptr);
#endif
#ifdef DEBUGMEMORY
inline void * chmalloc_memdebug(size_t size, char * filename, int lineNumber);
inline void * chrealloc_memdebug(void * oldptr, size_t size, char * filename, int lineNumber);
inline void * chcalloc_memdebug(size_t nmemb, size_t size, char * filename, int lineNumber);
inline void * chfree_memdebug(void *ptr, char * filename, int lineNumber);
#define chmalloc(size) chmalloc_memdebug(size, __FILE__, __LINE__)
#define chrealloc(oldptr, size) chrealloc_memdebug(oldptr, size, __FILE__, __LINE__)
#define chcalloc(nmemb, size) chcalloc_memdebug(nmemb, size, __FILE__, __LINE__)
#define chfree(ptr) chfree_memdebug(ptr, __FILE__, __LINE__)
#endif
inline FILE * chfopen(const char * path, const char * mode);
void timestamp(const char * msg);
void get_ts_msmc(double * ts, const int n);
void get_ts_psmc(double * ts, const double tmax, const double n);

#endif
