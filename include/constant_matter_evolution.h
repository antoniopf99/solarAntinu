#ifndef CONSTANT_MATTER_EVOLUTION_H
#define CONSTANT_MATTER_EVOLUTION_H

#include<iostream>
#include"putzer_algorithm.h"
#include"oscillation.h"
#include"earth_model.h"
#include"general_adiabatic.h"

class pauli_vector{
	public:
		complex<double> n[3];
		complex<double> dot_sigma[2][2];

		const complex<double> PAULI[3][2][2] = {
											   	{
											   		{0, 1},
											   		{1, 0}
											   	},
											   	{
											   		{0,-I},
											   		{I, 0}
											   	},
											   	{
											   		{1, 0},
											   		{0,-1}
											   	}
											   };
		void seta_dot_sigma();
};

double castle_wall_prob_2gen(double E,
	double L,
	double theta,
	double Dm,
	double T1,
	double T2,
	double n1,
	double n2);

double constant_matter_2gen(double E, double L, double theta, double Dm, double ne);

double constant_matter_prob_3gen(int nu_emitted, int nu_detected, double E, double L, double Ne, osc_param op);

double sun_earth_prob_2layers(double Rfrac, double E, osc_param op);

double sun_earth_2layers_Average8B(double E, osc_param op);

double sun_earth_prob_5layers(double Rfrac, double E, osc_param op);

double atmospheric_neutrino_5layers(int flavor_ini, int flavor_end, double E, osc_param op);

#endif