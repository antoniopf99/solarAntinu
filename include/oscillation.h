#ifndef OSCILATION_H
#define OSCILATION_H

#include<cmath>
#include<complex>
#include<fstream>
#include<gsl/gsl_errno.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_odeiv2.h>
#include"solarmodel.h"

using namespace std;

#ifndef IMAGINARY_UNIT
#define IMAGINARY_UNIT
const complex<double> I = complex<double>(0.0,1.0);
#endif //IMAGINARY_UNIT

/*
Class of neutrino state. Given a basis (e.g. vacuum, flavor), this class stores the probability amplitude
of finding the neutrino in each state of the base, identified with an index 0,1,2
*/
typedef struct nu_state{
	complex<double> index[3];
} nu_state;

/*debugging function*/
bool is_nan_nu_state(nu_state nu);

class osc_param{
	public:
		bool prob_ee;

		double theta12, theta13, theta23, delta_CP;
		double Dmq21, Dmq31;

		complex<double> U[3][3];

		double V[3][3];

		int n_exp_bins, i;
		double *exposure;
		double cosz;

		double h;

		void seta_osc_param(bool twogen);
		void constroi_PMNS();
};

/*initializes a PMNS matrix*/
void seta_PMNS(complex<double> U[][3], bool twogen);

/*
nu_em = EMITTED NEUTRINO, 
nu_ab = NEUTRINO THAT nu_em OSCILATTED TO,
x = DISTANCE TRAVELED IN km,
E = nu_em ENERGY/MOMENTUM IN GeV,
ordem_normal = IF YOU'RE USING NORMAL ORDERING OF THE MASSES.
*/
double survivalProb(int nu_em, int nu_ab, double x, double E, bool ordem_normal);

/*Now considering only one generation*/
double survivalProb1GEN(long double x, long double E, long double Dmq, long double theta, bool bestfit);

double survivalProbConstMatter(long double x, long double E, long double Dmq, long double theta12, long double theta13);

/*
nu_type_production = EMITTED NEUTRINO, 
nu_type_measure = NEUTRINO THAT nu_type_production OSCILATTED TO,
L = DISTANCE TRAVELED IN km,
E_nu = nu_em ENERGY/MOMENTUM IN GeV,
matter_density = MATTER DENSITY IN g/cm³
*/
double survival_prob_general(int nu_type_production, int nu_type_measure, double matter_density, long double E_nu, long double L);

/*
INCORRECT DESCRIPTION, UPDATE NEEDED
returns the survival probability of a solar neutrino as a function of neutrino energy and parameters
	 seta_dados_Ne() MUST BE CALLED
*/
double survival_prob_sun(double E_nu, double x0, void *param);

//seta_dados_Ne MUST BE CALLED
double survivalProb_2genSun_Ad(double E, double Rfrac);

//seta_dados_Ne MUST BE CALLED
double survivalProb_2genSun_Ad_Boron(double E);

//seta_dados_Ne MUST BE CALLED
double survivalProb_3genSun_Ad(double E, double Rfrac);

//seta_dados_Ne MUST BE CALLED
double survivalProb_3genSun_Ad_freeparam(double E, double Rfrac, osc_param op);

//seta_dados_Ne MUST BE CALLED
double survivalProb_3genSun_Ad_Boron(double E);

#endif //OSCILATION_H