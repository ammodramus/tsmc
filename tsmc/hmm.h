#ifndef HMM_H
#define HMM_H

#include "definitions.h"

typedef struct 
{
    int n; // number of changepoints
    int numStates;
    double * ts;
    double maxT;
    double * lambdas;
    double * deltas;
    double * intervalOmegas;
    double * pis;
    double * Eijs3s;
    double * Eijs2s;
    double ** qts;
    double ** pts;
    double rho;
    double theta;
    double Td;
    fourd * emissions;
} Hmm;

inline int get_index(int i, int j, int n);
void Hmm_init(Hmm * hmm, const int n);
void Hmm_free(Hmm * hmm);
void Hmm_set_lambdas(Hmm * hmm, const int n, const double * lambdas);
void Hmm_set_ts_and_deltas(Hmm * hmm, const double * ts);
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
void Hmm_get_pts(Hmm * hmm);
void Hmm_get_emissions(Hmm * hmm);
inline void Hmm_set_theta(Hmm * hmm, double theta);
inline void Hmm_set_rho(Hmm * hmm, double rho);
inline void Hmm_set_Td(Hmm * hmm, double Td);
void Hmm_make_hmm(Hmm * hmm, double * lambdas, double * ts,
    int numChangepoints, double theta, double rho, double Td, int * error);
void Hmm_print_pis(Hmm * hmm);
void Hmm_print_demography(Hmm * hmm);
void Hmm_print_expectations(Hmm * hmm);
void Hmm_print_qts(Hmm * hmm);
void Hmm_print_pts(Hmm * hmm);
void Hmm_print_emissions(Hmm * hmm);

#endif
