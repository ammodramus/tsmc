#ifndef HMM_H
#define HMM_H

typedef struct 
{
    int n; // number of changepoints
    int numStates;
    double * ts; // 
    double * lambdas;
    double * deltas;
    double * intervalOmegas;
    double * pis;
    double * Eijs3s;
    double * Eijs2s;
    double ** qts;
} Hmm;

inline int get_index(int i, int j, int n);
void Hmm_init(Hmm * hmm, int n, double * ts);
void Hmm_free(Hmm * hmm);
void Hmm_set_lambdas(Hmm * hmm, const int n, const double * lambdas);
void Hmm_make_omega_intervals(Hmm * hmm);
double get_omega_interval_interval(Hmm * hmm, int a, int b);
void Hmm_get_pis(Hmm * hmm);
void Hmm_get_expectations(Hmm * hmm);
double get_omega_Es3_interval(Hmm * hmm, const int b, const int Es3_i, 
        const int Es3_j);
double get_omega_Es2_interval(Hmm * hmm, const int b, const int Es2_i, 
        const int Es2_j);
double get_omega_interval_Es3(Hmm * hmm, const int a, const int Es3_i, 
        const int Es3_j);
double get_omega_interval_Es2(Hmm * hmm, const int a, const int Es2_i,
        const int Es2_j);
double get_omega_Es3_Es2(Hmm * hmm, const int Es3_i, const int Es3_j);
void Hmm_get_qts(Hmm * hmm);

#endif
