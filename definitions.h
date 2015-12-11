#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <stdlib.h>

#define DEBUGREPORTI(x) fprintf(stderr, #x " = %i\n", x)
#define DEBUGREPORTF(x) fprintf(stderr, #x " = %f\n", x)

typedef double threed[3];

void * chmalloc(size_t size);

#endif
