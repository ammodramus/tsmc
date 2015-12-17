#ifndef EH_H
#define EM_H

typedef struct 
{
    Data * dat;
    Hmm hmm[2];
    int hmmFlag;
    int numSeqs;
    int numHmmStates;
    int numIterations;
    int maxIterations;
    SeqType seqtype;
    double *** forward;
    double *** backward;
    double *** gamma;
    double ** expect;
    double ** normConst;
} Em;

void Em_init(Em * em, Data * dat, double * lambdas, double * ts, int n, 
        double initTheta, double initRho, double initTd, int maxIterations);
void Em_free(Em * em);
void Em_get_forward(Em * em);
void Em_get_backward(Em * em);
void Em_get_expectations(Em * em);
void Em_iterate(Em * em);
double Em_get_expected_log_likelihood(Em * em, const int hmmIdx);
double objective_function(double * par);
void Em_print_forward(Em * em);
void Em_print_backward(Em * em);
void Em_print_norm_const(Em * em);
void Em_print_gamma(Em * em);
void Em_print_expect(Em * em);

#endif
