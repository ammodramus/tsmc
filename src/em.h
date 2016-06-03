#ifndef EH_H
#define EM_H

#include "options.h"

typedef struct 
{
    Data * dat;
    Hmm hmm[2];
    int hmmFlag;
    int numSeqs;
    int flagDt;
    int numHmmStates;
    int curIteration;
    int maxIterations;
    int numOptimizations;
    int numFreeLambdas;
    SeqType seqtype;
    double *** forward;
    double *** backward;
    double *** gamma;
    double ** expectTransitions;
    sixd * expectEmissions;
    double ** normConst;
    double * freeLambdas;
    int * lambdaCounts;
    int asexEnabled;
} Em;

void Em_init(Em * em, Data * dat, double * ts,
        double initRho, int numFreeLambdas, 
        int n, int numEmIterations, int numOptimizations, int * lambdaCounts,
        int asexEnabled, int diptripflag);
void Em_free(Em * em);
double Em_get_initial_rho(Data * dat);
void Em_get_initial_theta_and_Td(Em * em, double * out);
void Em_get_forward(Em * em);
void Em_get_backward(Em * em);
void Em_get_expectations(Em * em);
void Em_iterate(Em * em);
void Em_iterate_asex(Em * em);
void Em_iterate_dt(Em * em);
void Em_iterate_no_asex(Em * em);
double Em_get_expected_log_likelihood(Em * em, const int hmmIdx);
double objective_function_asex(double * par);
double objective_function_no_asex(double * par);
double objective_function_dt(double * par);
void Em_print_forward(Em * em);
void Em_print_backward(Em * em);
void Em_print_norm_const(Em * em);
void Em_print_gamma(Em * em);
void Em_print_expect(Em * em);
double Em_get_loglikelihood(Em * em);
void Em_print_iteration(Em * em);

#endif
