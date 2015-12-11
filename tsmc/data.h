#ifndef DATA_H
#define DATA_H

#include "definitions.h"

/* Notes on data organization
 * ==========================
 *
 * 
 * Data struct needs:
 * -------------------
 * - number of seqs (int)
 * - type of data (polarized or not) (int or typedef'ed enum)
 * - data (Seq *)
 * - current number of allocated slots (int)
 *
 * Data methods/functions:
 * -----------------------
 *   - read data
 *   - allocate space for storing chroms
 *   - allocate a new Seq
 * 
 * Seq struct needs:
 * -----------------
 * - length (int)
 * - maxlength (int)
 * - index (int)
 * - data (char *)
 *
 */

typedef struct
{
    int len;
    int maxLen;
    int index;
    char * data;
} Seq;

typedef struct
{
    int numSeqs;
    int maxNumSeqs;
    int polarized;
    Seq * chroms;
} Data;

#endif
