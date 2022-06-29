#include<iostream>
#include"magnus_evolution.h"

nu_state magnus_evolution(nu_state nu_ini,
	double lenght_scale,
	double beginning_of_tragectory,
	double end_of_tragectory,
	double E,
	double (*matter_density)(double, void*),
	void* density_parameters,
	osc_param op)
{
	/*the vacuum portion H_0*/
	double H_0[3][3];
	double cte = 2533.8653588396974*lenght_scale; /* 1e-12 * 1/(2*hc) in units of 1/MeVkm */
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			H_0[i][j] = 0.0;
	H_0[1][1] = cte*op.Dmq21/E;
	H_0[2][2] = cte*op.Dmq31/E;

	/*the matter potential, in units of sqrt(2)*G_F*n_e*/
	complex<double> W[3][3];
	double V[3][3];
	double cte2 = 0.00038679283503942005*lenght_scale; /* 1e-6 sqrt(2)*GF */
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			W[i][j] = 0.0;
			V[i][j] = cte2*op.V[i][j];
		}
	}

	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			for(int k = 0; k < 3; k++) 
				W[i][j] += conj(op.U[k][i])*V[k][k]*op.U[k][j]; /*assuming that V is diagonal!*/

	/*************/
	/*COMMUTATORS*/
	/*************/
	complex<double> H_0cW[3][3];
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			H_0cW[i][j] = 0;
	
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			for(int k = 0; k < 3; k++)
				H_0cW[i][j] += H_0[i][k]*W[k][j] - W[i][k]*H_0[k][j];

	/*******************/
	/*INITIAL CONDITION*/
	/*******************/
	complex<double> Vatual[3], Vprox[3];
	for(int i = 0; i < 3; i++)
		Vatual[i] = nu_ini.index[i];

	/****************************/
	/*PERFORMING THE INTEGRATION*/
	/****************************/
	complex<double> Omega[3][3];
	double xi_plus, xi_minus, v_plus, v_minus, h = 0.00005, sqrt3 = sqrt(3);

	for(double xi = beginning_of_tragectory; xi < end_of_tragectory; xi += h){
		for(int i = 0; i < 3; i++)
			for(int j = 0; j < 3; j++)
				Omega[i][j] = 0.0;
		
		/****************************/
		/*COMPUTING THE OMEGA MATRIX*/
		/****************************/
		xi_plus = xi + (1.0 + (1.0/sqrt3))*(h/2.0);
		xi_minus = xi + (1.0 - (1.0/sqrt3))*(h/2.0);

		v_plus = matter_density(xi_plus, density_parameters);
		v_minus = matter_density(xi_minus, density_parameters);

		/*i include an extra -i factor to make it hermitian*/
		for(int i = 0; i < 3; i++){
			for(int j = 0; j < 3; j++){
				Omega[i][j] = -I*(-I*(H_0[i][j] + 0.5*(v_plus+v_minus)*W[i][j])*h + (sqrt3/12.0)*(v_plus-v_minus)*H_0cW[i][j]*h*h);
			}
		}

		/********************************************/
		/*MATRIX EXPONENTIAL WITH PUTZER'S ALGORITHM*/
		/********************************************/
		matrixexp_times_vector_3x3_putzer(Omega, Vatual, Vprox);

		/*setting back current step*/
		for(int i = 0; i < 3; i++)
			Vatual[i] = Vprox[i];
	}

	nu_state nu;
	for(int i = 0; i < 3; i++)
		nu.index[i] = Vatual[i];

	return(nu);
}

double sun_Ne(double Rfrac, void *p){
	return(pow(10.0, electron_density_sun_NA(Rfrac)));
}

double survival_prob_sun_magnus(double Rfrac, double E, osc_param op){
	/***************************************/
	/*EVALUATING THE AVERAGED SURVIVAL PROB*/
	/***************************************/

	/*Sun's radius*/
	double Rs = 696340; /*in km*/

	/*Initial condition*/
	nu_state nu_ini;
	for(int i = 0; i < 3; i++)
		nu_ini.index[i] = op.U[0][i]; /*0, electron neutrino!*/

	/*Evolution with the Magnus expansion*/
	void *density_parameters = NULL;
	nu_state nu = magnus_evolution(nu_ini, Rs, Rfrac, 1, E, &sun_Ne, density_parameters, op);

	/*Incoherent average survival probability*/
	double resultado = 0;
	for(int i = 0; i < 3; i++)
		resultado += norm(op.U[0][i])*norm(nu.index[i]);

	return(resultado);
}

double survival_prob_sun_magnus_Average8B(double E, osc_param op){
	/*points used to compose the average of boron neutrino production points*/
	int N = 25;
	double resultado = 0, normalizacao = 0;
	double RfracMax = 0.2;
	double h = RfracMax/N;
	double boronFrac;
	for(double Rfrac = 0; Rfrac < RfracMax; Rfrac += h){
		boronFrac = interpola_dados_solares(Rfrac, 8);
		resultado += boronFrac*survival_prob_sun_magnus(Rfrac, E, op);
		normalizacao += boronFrac;
	}
	return(resultado/normalizacao);
}

double earth_Ne(double x, void* p){
	double* cosz = (double*)p;
	return(earth_electron_density_NA(radial_distance_core(*cosz, 0.0, x)));
}

double survival_prob_sun_earth_magnus(double Rfrac, double E, osc_param op){
	if(op.cosz <= 0)
		return(survival_prob_sun_magnus(Rfrac, E, op));

	/*Sun's and earth's radius*/
	double Rs = 696340; /*in km*/
	double Re = 6371; /*in km*/

	double* ptr_cosz = &(op.cosz);
	void* density_parameters = (void*)ptr_cosz;

	nu_state nu_ini, nu_decoh, nu_earth_evolved;
	complex<double> amp_e;
	double Pe[3] = {0, 0, 0};

	/*initial condition for sun evolution*/
	for(int i = 0; i < 3; i++)
		nu_ini.index[i] = op.U[0][i]; /*0, electron neutrino!*/

	/*evolving inside the sun*/
	nu_decoh = magnus_evolution(nu_ini, Rs, Rfrac, 1, E, &sun_Ne, density_parameters, op);

	/*now we evolve the 1 and 2 matter states inside the earth (and use unitarity for the third)*/
	for(int i = 0; i < 2; i++){

		/*initial condition for earth evolution, matter eigenstate*/
		for(int j = 0; j < 3; j++){
			if(j == i)
				nu_ini.index[j] = 1.0;
			else
				nu_ini.index[j] = 0.0;
		}

		/*earth evolution*/
		nu_earth_evolved = magnus_evolution(nu_ini, Re, 0, distance_inside_earth(op.cosz, op.h), E, &earth_Ne, density_parameters, op);
			
		/*survival probability*/
		amp_e = 0.0;
		for(int j = 0; j < 3; j++)
			amp_e += op.U[0][j]*nu_earth_evolved.index[j];
		Pe[i] = norm(amp_e);
	}

	Pe[2] = (1 - Pe[0] - Pe[1]);

	/*the survival probs with relative admixtures*/
	return(norm(nu_decoh.index[0])*Pe[0] + norm(nu_decoh.index[1])*Pe[1] + norm(nu_decoh.index[2])*Pe[2]);
}

double survival_prob_sun_earth_magnus_Average8B(double E, osc_param op){
	/*points used to compose the average of boron neutrino production points*/
	int N = 20;
	double resultado = 0, normalizacao = 0;
	double RfracMax = 0.2;
	double h = RfracMax/N;
	double boronFrac;
	for(double Rfrac = 0; Rfrac < RfracMax; Rfrac += h){
		boronFrac = interpola_dados_solares(Rfrac, 8);
		resultado += boronFrac*survival_prob_sun_earth_magnus(Rfrac, E, op);
		normalizacao += boronFrac;
	}
	return(resultado/normalizacao);
}

double survival_prob_sun_ad_earth_magnus(double Rfrac, double E, osc_param op){
	if(op.cosz <= 0)
		return(survival_prob_sun_general_adiabatic(Rfrac, E, op));

	//earth's radius
	double Re = 6371; //in km

	double* ptr_cosz = &(op.cosz);
	void* density_parameters = (void*)ptr_cosz;

	nu_state nu_ini, nu_decoh, nu_earth_evolved;
	complex<double> amp_e;
	double Pe[3] = {0, 0, 0}, eigenv[3];

	//adiabatic decoherent evolution inside the sun
	nu_decoh = sun_adiabatic_decoherent_evolution(Rfrac, E, op, eigenv);

	//now we evolve the 1 and 2 matter states inside the earth (and use unitarity for the third)
	for(int i = 0; i < 2; i++){

		//initial condition for earth evolution, matter eigenstate
		for(int j = 0; j < 3; j++){
			if(j == i)
				nu_ini.index[j] = 1.0;
			else
				nu_ini.index[j] = 0.0;
		}

		//earth evolution
		nu_earth_evolved = magnus_evolution(nu_ini, Re, 0, distance_inside_earth(op.cosz, op.h), E, &earth_Ne, density_parameters, op);
			
		//survival probability
		amp_e = 0.0;
		for(int j = 0; j < 3; j++)
			amp_e += op.U[0][j]*nu_earth_evolved.index[j];
		Pe[i] = norm(amp_e);
	}

	Pe[2] = (1 - Pe[0] - Pe[1]);

	//the survival probs with relative admixtures
	return(norm(nu_decoh.index[0])*Pe[0] + norm(nu_decoh.index[1])*Pe[1] + norm(nu_decoh.index[2])*Pe[2]);
}

double survival_prob_sun_ad_earth_magnus_Average8B(double E, osc_param op){
	//points used to compose the average of boron neutrino production points
	int N = 20;
	double resultado = 0, normalizacao = 0;
	double RfracMax = 0.2;
	double h = RfracMax/N;
	double boronFrac;
	for(double Rfrac = 0; Rfrac < RfracMax; Rfrac += h){
		boronFrac = interpola_dados_solares(Rfrac, 8);
		resultado += boronFrac*survival_prob_sun_ad_earth_magnus(Rfrac, E, op);
		normalizacao += boronFrac;
	}
	return(resultado/normalizacao);
}