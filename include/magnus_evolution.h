#ifndef MAGNUS_EVOLUTION_H
#define MAGNUS_EVOLUTION_H

#include<cmath>
#include<complex>
#include"putzer_algorithm.h"
#include"solarmodel.h"
#include"oscillation.h"
#include"earth_model.h"
#include"general_adiabatic.h"

using namespace std;

/*
This funtion returns the desired evolved neutrino state
*/
nu_state magnus_evolution(nu_state nu_ini,
	double lenght_scale,
	double beginning_of_tragectory,
	double end_of_tragectory,
	double E,
	double (*matter_density)(double, void *),
	void *density_parameters,
	osc_param op);
/*
Returns the survival probability of a solar electron neutrino produced at Rfrac, in units
of solar radius.
E in MeV
V_e, V_mu and V_tau in units of sqrt(2)*G_F*n_e
*/
double survival_prob_sun_magnus(double Rfrac, double E, osc_param op);

//seta_dados_Ne MUST BE CALLED
double survival_prob_sun_magnus_Average8B(double E, osc_param op);

double survival_prob_sun_earth_magnus(double Rfrac, double E, osc_param op);

double survival_prob_sun_earth_magnus_Average8B(double E, osc_param op);

double survival_prob_sun_ad_earth_magnus(double Rfrac, double E, osc_param op);

#endif
