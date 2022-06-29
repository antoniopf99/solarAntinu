#ifndef STATISTICS_H
#define STATISTICS_H

#include<cmath>
#include<gsl/gsl_multimin.h>

typedef struct STAT_PARAMS{
	int nBins;
	double *N_obs;
	double *T;
	double *BG;
	double sd;
	double sb;
} stat_params;

double chisquare_twopulls_notminimized(const gsl_vector *x, void *params);

double chisquare_twopulls_JUNO(stat_params p, bool verbo);

/*this function returns a simple gaussian chi square, using the square root of
N_obs as the incertanty (as if the data were poissonically distributed)*/
double simple_chisquare(int nBins, double N_pred[], double N_obs[]);

#endif