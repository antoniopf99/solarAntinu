#ifndef DETECTORS_H
#define DETECTORS_H

#include<iostream>
#include<cmath>
#include<cstring>
#include<fstream>
#include<vector>
#include<gsl/gsl_sf_erf.h>
#include<gsl/gsl_integration.h>
#include"crosssections.h"
#include"oscillation.h"

using namespace std;

/*The absolute and relative errors used in GSL's numeric integration*/
const double ERRO_ABS = 1e-1;//1e-1
const double ERRO_REL = 1e-4;//1e-4

/*The maximum number of subintervals used in some GSL integration routines*/
const int MAX_SUB_INT = 500;

typedef struct SIGNAL_PARAM{
	double (*flux)(double);
	double max_flux_E;
	double (*xsection)(double, double, void *);
	void *xsection_param;
	double (*oscillation_prob)(double, void *);
	void *oscillation_param;
} signal_param;

typedef struct DETECTOR_PARAM{
	int numberOfBins;
	double Ebinmin;
	double Ebinmax;
	double binsize;
	double Eres;
	double (*Ntau)(double);
} detector_param;

typedef struct INT_PARAM{
	signal_param sig;
	detector_param detec;
	double Tanalmin;
	double Tanalmax;
	int N; /*number of interpolation points*/
	double *T; /*recoil energy points used for interpolation*/
	double *countRate; /*count rate points used for interpolation*/
	double E_current_bin; /*minimum energy of the current bin*/
} integration_parameters;

double events_in_bin(signal_param sig, detector_param detec);

double countRate_nu_e_Scattering_integrand(double E, void *sig_param);

double countRate_nu_e_Scattering(double T, signal_param sig);

double calcula_eventos_Bins_integrand(double T, void *quad_param);

void calcula_eventos_Bins(double *Bins, signal_param sig, detector_param detec);

void seta_BG(double *BG, int numberOfBins, string bg_file_name);

void seta_exposure(double *exposure, int numberOfExpBins, string exp_file_name);

#endif