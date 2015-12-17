#ifndef OPTIONS_H
#define OPTIONS_H

#define MAXFILENAMESIZE 200
#define MAXNUMFREELAMBDAS 200

typedef struct
{
    int n;
    char filename[MAXFILENAMESIZE];
    int numFreeLambdas;
    int lambdaCounts[MAXNUMFREELAMBDAS];
    int numEmIterations;
    char paramString[200];
} Options;

#endif
