#ifndef GENERAL_ADIABATIC_H
#define GENERAL_ADIABATIC_H

#include<gsl/gsl_matrix_double.h>
#include<gsl/gsl_eigen.h>
#include"oscillation.h"
#include"solarmodel.h"

nu_state sun_adiabatic_decoherent_evolution(double Rfrac, double E, osc_param op, double eigenv[3]);

double survival_prob_sun_general_adiabatic(double Rfrac, double E, osc_param op);

double masses_values(int index, double Rfrac, double E, double V_e, double V_mu, double V_tau);

double survival_prob_sun_adiabatic_Zplike(double Rfrac, double E, double V_e, double V_mu, double V_tau);

#endif