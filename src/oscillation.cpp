#include<iostream>
#include"oscillation.h"

bool is_nan_nu_state(nu_state nu){
	bool ehnan = false;
	for(int i = 0; i < 3 && !ehnan; i++)
		if(isnan(real(nu.index[i])) || isnan(imag(nu.index[i])))
			ehnan = true;
	return(ehnan);
}

void zera_H_c(complex<double> M[][3]){
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			M[i][j] = 0;
}

void osc_param::seta_osc_param(bool twogen){
	theta12 = asin(sqrt(0.307));//34.5 * M_PI / 180;
	if(twogen){
		theta13 = 0;
		theta23 = 0;
	}
	else{
		theta13 = asin(sqrt(0.0215));//8.45 * M_PI / 180;
		theta23 = asin(sqrt(0.425));//47.7 * M_PI / 180;
	}
	delta_CP = 1.38*M_PI;//195.0 * M_PI / 180; /*195*/

	Dmq21 = 4.8;//7.55e-5;
	Dmq31 = 2.56;//2.5e-3;

	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			V[i][j] = 0;

	V[0][0] = 1;
}

void osc_param::constroi_PMNS(){
	double s12 = sin(theta12);
	double s13 = sin(theta13);
	double s23 = sin(theta23);
	double c12 = cos(theta12);
	double c13 = cos(theta13);
	double c23 = cos(theta23);

	complex<double> faseCP = exp(I * delta_CP);

	U[0][0] = c12 * c13;
	U[0][1] = s12 * c13;
	U[0][2] = s13 * conj(faseCP);
	U[1][0] = -s12 * c23 - c12 * s23 * s13 * faseCP;
	U[1][1] = c12 * c23 - s12 * s23 * s13 * faseCP;
	U[1][2] = s23 * c13;
	U[2][0] = s12 * s23 - c12 * c23 * s13 * faseCP;
	U[2][1] = -c12 * s23 - s12 * c23 * s13 * faseCP;
	U[2][2] = c23 * c13;
}

/******************/
/*PDG CONVENTION!!*/
/******************/
void seta_PMNS(complex<double> U[][3], bool twogen){

	/************************/
	/*PMNS MATRIX PARAMETERS*/
	/************************/
	double theta12 = 33.45 * M_PI / 180;
	double theta13;
	double theta23;
	if(twogen){
		theta13 = 0;
		theta23 = 0;
	}
	else{
		theta13 = 8.62 * M_PI / 180;
		theta23 = 42.1 * M_PI / 180;
	}

	double deltaCP = 195.0 * M_PI / 180; /*195*/

	double s12 = sin(theta12);
	double s13 = sin(theta13);
	double s23 = sin(theta23);
	double c12 = cos(theta12);
	double c13 = cos(theta13);
	double c23 = cos(theta23);

	complex<double> faseCP = exp(I * deltaCP);

	U[0][0] = c12 * c13;
	U[0][1] = s12 * c13;
	U[0][2] = s13 * conj(faseCP);
	U[1][0] = -s12 * c23 - c12 * s23 * s13 * faseCP;
	U[1][1] = c12 * c23 - s12 * s23 * s13 * faseCP;
	U[1][2] = s23 * c13;
	U[2][0] = s12 * s23 - c12 * c23 * s13 * faseCP;
	U[2][1] = -c12 * s23 - s12 * c23 * s13 * faseCP;
	U[2][2] = c23 * c13;
}

double survivalProb(int nu_em, int nu_ab, double x, double E, bool ordem_normal){

	complex<double> U[3][3];
	seta_PMNS(U,false);

	/*********************************************/
	/*NEUTRINO OSCILATION PARAMETERS, m^2 in eV^2*/
	/*********************************************/
	double Dmq21 = 7.42 * pow(10,-5);/*7.42 * pow(10,-5);*/

	double mq[3];

	/*Only differences between the masses are relevant, let's set the lightest to 0*/
	if(ordem_normal){
		/***********************/
		/*Dmq31 NORMAL ORDERING*/
		/***********************/
		double Dmq31 = 2.514 * pow(10,-3);/*2.514 * pow(10,-3);*/ /*Dmq31, m1 < m2 < m3*/

		mq[0] = 0; /*in the normal ordering, m1 is the lightest*/
		mq[1] = Dmq21; /*because mq1 = 0*/
		mq[2] = Dmq31;
	}
	else{
		/************************/
		/*Dmq32 INVERSE ORDERING*/
		/************************/
		double Dmq32 = -2.497 * pow(10,-3); /*Dmq32, m3 < m1 < m2*/
		
		mq[2] = 0; /*in the inverse ordering, m3 is the lightest*/
		mq[1] = - Dmq32;
		mq[0] = - Dmq21 - Dmq32;
	}

	complex<double> soma = 0;

	double conts = 2.533865358274169; /*1 / hbar (eVs) * c (m/s) * G/k*/

	/*So xoverE is in km/GeV*/
	for(int i = 0; i < 3; i++)
		soma += conj(U[nu_ab][i]) * U[nu_em][i] * exp(-I * conts * mq[i] * (x / E));
	return(norm(soma));
}

/*this Hamiltonian is in the mass basis!!*/
int Hamiltoniano(double x, const double nu_m[], double f[], void *parametros){

	/*recovering the parameters*/
	double *param = (double*)parametros;
	double E = param[0];
	double matter_density = param[1];

	/*mass mixing parameters*/
	double Dmq21 = 7.42e-5;
	double Dmq31 = 2.514e-3;
	
	/*************/
	/*HAMILTONIAN*/
	/*************/
	double H[6][6];
	complex<double> preH[3][3];
	complex<double> matterH[3][3];
	complex<double> U[3][3];
	complex<double> Udagger[3][3];
	seta_PMNS(U, false);
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			Udagger[i][j] = conj(U[j][i]);
	zera_H_c(preH);
	zera_H_c(matterH);

	double cte = 1.26693276;
	cte *= 2;
	double cte2 = 0.00019257376961400582;
	
	/*hamiltonian contribution from neutrino masses in mass basis*/
	preH[1][1] += (cte/E)*Dmq21;
	preH[2][2] += (cte/E)*Dmq31;

	/*matter potential hamiltonian in the mass basis*/
	matterH[0][0] += cte2*matter_density;

	/*adding contributions from matter potential in mass basis*/
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			for(int k = 0; k < 3; k++)
				for(int l = 0; l < 3; l++)
					preH[i][j] += Udagger[i][k]*matterH[k][l]*U[l][j];


	/*setting the six dimensional rela hamiltonian*/
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			H[i][j+3] = real(preH[i][j]);
			H[i+3][j] = -real(preH[i][j]);
			H[i][j] = imag(preH[i][j]);
			H[i+3][j+3] = imag(preH[i][j]);
		}
	}

	for(int i = 0; i < 6; i++){
		f[i] = 0;
		for(int j = 0; j < 6; j++){
			f[i] += (H[i][j]*nu_m[j]);
		}
	}
	
	return GSL_SUCCESS;
}

/*this Hamiltonian is the flavor basis!!*/
int Hamiltoniano_sol(double x, const double nu_m[], double f[], void *parametros){

	/*recovering the parameters*/
	double *param = (double*)parametros;
	double E = param[0];

	/*mass mixing parameters*/
	double Dmq21 = 7.42e-5;
	double Dmq31 = 2.514e-3;
	
	/*************/
	/*HAMILTONIAN*/
	/*************/
	double H[6][6];
	complex<double> preH[3][3];
	complex<double> matterH[3][3];
	complex<double> U[3][3];
	complex<double> Udagger[3][3];
	seta_PMNS(U, false);
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			Udagger[i][j] = conj(U[j][i]);
	zera_H_c(preH);
	zera_H_c(matterH);

	double cte = 1.26693276;
	cte *= 2;
	double cte2 = 3.8679283503942017e-07;
	
	/*hamiltonian contribution from neutrino masses in mass basis*/
	preH[1][1] += (cte/E)*Dmq21;
	preH[2][2] += (cte/E)*Dmq31;

	/*matter potential hamiltonian in the mass basis*/
	matterH[0][0] += cte2*pow(10.0 ,electron_density_sun_NA(x));

	/*adding contributions from matter potential in mass basis*/
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			for(int k = 0; k < 3; k++)
				for(int l = 0; l < 3; l++)
					preH[i][j] += Udagger[i][k]*matterH[k][l]*U[l][j];


	/*setting the six dimensional rela hamiltonian*/
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			H[i][j+3] = real(preH[i][j]);
			H[i+3][j] = -real(preH[i][j]);
			H[i][j] = imag(preH[i][j]);
			H[i+3][j+3] = imag(preH[i][j]);
		}
	}

	for(int i = 0; i < 6; i++){
		f[i] = 0;
		for(int j = 0; j < 6; j++){
			f[i] += (H[i][j]*nu_m[j]);
		}
	}
	
	return GSL_SUCCESS;
}

double survival_prob_general(int nu_type_production, int nu_type_measure, double matter_density, long double E_nu, long double L){
	
	/*numeric integration parameters*/
	const int dim = 6;
	const double erro_abs = 1e-8;
	const double erro_rel = 1e-8;
	int status;

	/*parameters of the evolution equation*/
	double parametros[2] = {(double)E_nu, matter_density};

	/*variables used in the numeric integration*/
	double h = 1e-6; /*passo numérico inicial*/


	/****************************/
	/*SETTING INITIAL CONDITIONS*/
	/****************************/
	/*defining the PMNS matrix and its dagger*/
	complex<double> U[3][3];
	complex<double> Udagger[3][3];
	seta_PMNS(U,false);
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			Udagger[i][j] = conj(U[j][i]);

	/*state vector in mass basis, first (last) 3 entries are the real (imaginary) parts*/
	double nu_m[6] = {0, 0, 0, 0, 0, 0};

	/*state vector in flavor space, complex*/
	complex<double> nu_f[3] = {0, 0, 0};

	/*initial condition*/
	nu_f[nu_type_production] = 1;

	/*converting initial condition to the mass basis 6 component vector*/
	for(int i = 0; i < 3; i ++){
		for(int j = 0; j < 3; j++){
			nu_m[i] += real(Udagger[i][j]*nu_f[j]);
			nu_m[3+i] += imag(Udagger[i][j]*nu_f[j]);
		} 
	}
	/****************************/

	/*integration interval*/
	double x_ini = 0;
	double x_fim = L;

	/*setting up the GSL machinary*/
	gsl_odeiv2_system sistema;

	sistema.function = Hamiltoniano;
	sistema.dimension = dim;
	sistema.params = parametros;

	gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(&sistema, gsl_odeiv2_step_rkf45, h, erro_abs, erro_rel);

	/*numeric integration*/
	status = gsl_odeiv2_driver_apply(driver, &x_ini, x_fim, nu_m);

	if(status != GSL_SUCCESS)
		cout << "ERRO!!!" << endl;

	/*converting back to the flavor basis*/
	complex<double> nu_m_complex[3];
	for(int i = 0; i < 3; i++)
		nu_m_complex[i] = nu_m[i] + I*nu_m[3+i];

	for(int i = 0; i < 3; i++){
		nu_f[i] = 0;
		for(int j = 0; j < 3; j++)
			nu_f[i] += U[i][j]*nu_m_complex[j];
	}

	double resultado = norm(nu_f[nu_type_measure]);
	gsl_odeiv2_driver_free (driver);
	return(resultado);
}

double survival_prob_sun(double E_nu, double x0, void *param){
		
	/*numeric integration parameters*/
	const int dim = 6;
	const double erro_abs = 1e-8;
	const double erro_rel = 1e-8;
	int status;

	/*parameters of the evolution equation*/
	double parametros[1] = {(double)E_nu};

	/*variables used in the numeric integration*/
	double h = 1; /*passo numérico inicial*/


	/****************************/
	/*SETTING INITIAL CONDITIONS*/
	/****************************/
	/*defining the PMNS matrix and its dagger*/
	complex<double> U[3][3];
	complex<double> Udagger[3][3];
	seta_PMNS(U,false);
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			Udagger[i][j] = conj(U[j][i]);

	/*state vector in mass basis, first (last) 3 entries are the real (imaginary) parts*/
	double nu_m[6] = {0, 0, 0, 0, 0, 0};

	/*state vector in flavor space, complex*/
	complex<double> nu_f[3] = {0, 0, 0};

	/*initial condition*/
	nu_f[0] = 1;

	/*converting initial condition to the mass basis 6 component vector*/
	for(int i = 0; i < 3; i ++){
		for(int j = 0; j < 3; j++){
			nu_m[i] += real(Udagger[i][j]*nu_f[j]);
			nu_m[3+i] += imag(Udagger[i][j]*nu_f[j]);
		} 
	}
	/****************************/

	/*integration interval*/
	double x_ini = x0;
	double x_fim = 696340000;

	/*setting up the GSL machinary*/
	gsl_odeiv2_system sistema;

	sistema.function = Hamiltoniano;
	sistema.dimension = dim;
	sistema.params = parametros;

	gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(&sistema, gsl_odeiv2_step_rkf45, h, erro_abs, erro_rel);

	/*numeric integration*/
	status = gsl_odeiv2_driver_apply(driver, &x_ini, x_fim, nu_m);

	if(status != GSL_SUCCESS)
		cout << "ERRO!!!" << endl;

	complex<double> nu_m_complex[3];
	for(int i = 0; i < 3; i++)
		nu_m_complex[i] = nu_m[i] + I*nu_m[3+i];

	double resultado = 0;
	for(int i = 0; i < 3; i++)
		resultado += norm(U[0][i])*norm(nu_m_complex[i]);
	gsl_odeiv2_driver_free(driver);
	return(resultado);
}

double survivalProb_2genSun_Ad(double E, double Rfrac){

	/*mixing angles*/
	double theta = 33.44 * M_PI / 180; /*theta12*/
	double si2 = sin(2*theta);
	double co2 = cos(2*theta);

	/*mass mixing parameter*/
	double Dmq = 7.42e-5; /*Dmq21*/

	double Ne = pow(10,electron_density_sun_NA(Rfrac));
	double A = 1.52649*pow(10,-13)*Ne*pow(10,6)*E;
	double DmqM = sqrt(pow(Dmq*co2 - A,2) + pow(Dmq*si2,2));
	double co2M = (Dmq*co2 - A)/(DmqM);

	return(0.5*(1 + co2M*co2));
}

double survivalProb_2genSun_Ad_Boron(double E){
	/*points used to compose the average of boron neutrino production points*/
	int N = 25;
	double resultado = 0, normalizacao = 0;
	double RfracMax = 0.2;
	double h = RfracMax/N;
	double boronFrac;
	for(double Rfrac = 0; Rfrac < RfracMax; Rfrac += h){
		boronFrac = interpola_dados_solares(Rfrac, 8);
		resultado += boronFrac*survivalProb_2genSun_Ad(E, Rfrac);
		normalizacao += boronFrac;
	}
	return(resultado/normalizacao);
}

double survivalProb_3genSun_Ad(double E, double Rfrac){
	/*mixing angles*/
	double theta = 33.45 * M_PI / 180; /*theta12*/
	double theta13 = 8.62 * M_PI / 180; /*theta13*/
	double si2 = sin(2*theta);
	double co2 = cos(2*theta);

	/*mass mixing parameter*/
	double Dmq = 7.42e-5; /*Dmq21*/

	double Ne = pow(10,electron_density_sun_NA(Rfrac));
	double A = 1.52649*pow(cos(theta13),2)*pow(10,-13)*Ne*pow(10,6)*E;
	double DmqM = sqrt(pow(Dmq*co2 - A,2) + pow(Dmq*si2,2));
	double co2M = (Dmq*co2 - A)/(DmqM);

	double P12 = 0.5*(1 + co2M*co2);

	return(pow(cos(theta13),4)*P12 + pow(sin(theta13),4));
}

double survivalProb_3genSun_Ad_freeparam(double E, double Rfrac, osc_param op){
	
	double si2 = sin(2*(op.theta12));
	double co2 = cos(2*(op.theta12));
	double Dmq = op.Dmq21;
	double theta13 = op.theta13;

	double Ne = pow(10,electron_density_sun_NA(Rfrac));
	double A = 1.52649*pow(cos(theta13),2)*pow(10,-13)*Ne*pow(10,6)*E;
	double DmqM = sqrt(pow(Dmq*co2 - A,2) + pow(Dmq*si2,2));
	double co2M = (Dmq*co2 - A)/(DmqM);

	double P12 = 0.5*(1 + co2M*co2);

	return(pow(cos(theta13),4)*P12 + pow(sin(theta13),4));
}

double survivalProb_3genSun_Ad_Boron(double E){
	/*points used to compose the average of boron neutrino production points*/
	int N = 15;
	double resultado = 0, normalizacao = 0;
	double RfracMax = 0.2;
	double h = RfracMax/N;
	double boronFrac;
	for(double Rfrac = 0; Rfrac < RfracMax; Rfrac += h){
		boronFrac = interpola_dados_solares(Rfrac, 8);
		resultado += boronFrac*survivalProb_3genSun_Ad(E, Rfrac);
		normalizacao += boronFrac;
	}
	return(resultado/normalizacao);
}

double survivalProb1GEN(long double x, long double E, long double theta, long double Dmq, bool bestfit){
	if(bestfit){

		/************************/
		/*PMNS MATRIX PARAMETERS*/
		/************************/
		long double theta12 = 0.642432; /*33.44 * M_PI / 180;*/

		/*********************************************/
		/*NEUTRINO OSCILATION PARAMETERS, Dm^2 in eV^2*/
		/*********************************************/
		long double Dmq21 = 7.42 * pow(10,-5);/*7.42 * pow(10,-5);*/

		return(1 - pow(sin(2*theta12),2) * pow(sin(1.2669 * ((Dmq21*x)/E)),2));
	}
	return(1 - pow(sin(2*theta),2) * pow(sin(1.2669 * ((Dmq*x)/E)),2));
}

double survivalProbConstMatter(long double x, long double E, long double Dmq, long double theta12, long double theta13){
	/*******************/
	/*MATTER PARAMETERS*/
	/*******************/
	/*Matter density in g/cm³*/
	double rho = 2.7;
	/*Relative number density*/
	double Ye = 0.5;
	
	/*E comes in GeV = 1e9eV, so there's a factor of 1e(-14+9) = 1e-5*/
	long double A = 2 * E * 7.6 * pow(10,-5) * pow(cos(theta13),2) * Ye * rho / Dmq;

	long double factor = pow(cos(2*theta12) + A, 2) + pow(sin(2*theta12), 2);

	long double sin2theta12M = pow(sin(2*theta12), 2) / factor;
	long double DmqM = Dmq * sqrt(factor);

	long double Peetilde = 1 - sin2theta12M * pow(sin(1.2669 * ((DmqM*x)/E)),2);
	
	return(pow(cos(theta13),4) * Peetilde + pow(sin(theta13),4));
}