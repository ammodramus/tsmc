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
    int numOptimizations;
    int asexEnabled;
    int flagDt;
    int psmcIntervals;
    double maxT;
    char paramString[200];
} Options;

void Options_parse_options(Options * opt, int argc, char ** argv);
void Options_print_help_statement();
void Options_set_defaults(Options * opt);
void Options_get_lambda_params(Options * opt);
void Options_parse_options(Options * opt, int argc, char ** argv);
void Options_print_options(Options * opt);

#endif
