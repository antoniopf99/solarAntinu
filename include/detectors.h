#ifndef DETECTORS_H
#define DETECTORS_H

#include<iostream>
#include<cmath>
#include<vector>
#include<gsl/gsl_sf_erf.h>

using namespace std;

double events_in_bin(
	long double (*Flux)(long double), 
	long double (*Xsection)(long double, long double, double *), 
	double *XsectionParameters, 
	long double (*survivalProb)(long double, double *), 
	double *survivalProbParameters,
	double Ebinmin,
	double binsize,
	long double Ntau,
	double Eres);

long double countRate_nu_e_Scattering(
	long double (*Flux)(long double), 
	long double (*Xsection)(long double, long double, double *), 
	double *XsectionParameters, 
	long double (*survivalProb)(long double, double *), 
	double *survivalProbParameters,
	double T);

void calcula_eventos_Bins(
	double *Bins,
	long double Ebinmin,
	long double Ebinmax,
	int numberOfBins,
	long double (*Ntau)(long double),
	long double (*Flux)(long double), 
	long double (*Xsection)(long double, long double, double *), 
	double *XsectionParameters, 
	long double (*survivalProb)(long double, double *), 
	double *survivalProbParameters,
	double Eres);

#endif