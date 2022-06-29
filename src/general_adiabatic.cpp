#include<iostream>
#include"general_adiabatic.h"

void printa_matriz(complex<double> M[][3]){
	cout << endl;
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++)
			cout << M[i][j] << " ";
		cout << endl;
	}
	cout << endl;
}

nu_state sun_adiabatic_decoherent_evolution(double Rfrac, double E, osc_param op, double eigenv[3]){
	/***************************************************/
	/*LET'S CONSTRUCT H IN UNITS OF (micro 1e-6) \mu eV*/
	/***************************************************/

	/*the vacuum portion H_0*/
	double H_0_matter[3][3];
	complex<double> H_0_flavor[3][3];
	double cte;
	cte = 1;
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			H_0_flavor[i][j] = H_0_matter[i][j] = 0.0;
	H_0_matter[1][1] = cte*op.Dmq21/(2*E);
	H_0_matter[2][2] = cte*op.Dmq31/(2*E);

	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			for(int k = 0; k < 3; k++)
				H_0_flavor[i][j] += op.U[i][k]*H_0_matter[k][k]*conj(op.U[j][k]);

	/*the matter potential, in units of sqrt(2)*G_F*n_e*/
	double Ne = pow(10.0, electron_density_sun_NA(Rfrac));
	double cte2, V[3][3];
	cte2 = 7.6324662e-8; /* 1e-6 sqrt(2)*GF */
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			V[i][j] = cte2*Ne*op.V[i][j];

	/*the full Hamiltonian*/
	complex<double> H[3][3];
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			H[i][j] = H_0_flavor[i][j] + V[i][j];

	/***************/
	/*GSL MACHINARY*/	
	/***************/
	gsl_eigen_hermv_workspace *workspace = gsl_eigen_hermv_alloc(3);
	gsl_matrix_complex *gsl_H = gsl_matrix_complex_alloc(3, 3);
	gsl_matrix_complex *gsl_U_matter = gsl_matrix_complex_alloc(3, 3);
	gsl_vector *gsl_eigenvalues = gsl_vector_alloc(3);

	gsl_complex z;
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			GSL_SET_COMPLEX(&z, real(H[i][j]), imag(H[i][j]));
			gsl_matrix_complex_set(gsl_H, i, j, z);
		}
	}

	gsl_eigen_hermv(gsl_H, gsl_eigenvalues, gsl_U_matter, workspace);
	gsl_eigen_hermv_sort(gsl_eigenvalues, gsl_U_matter, GSL_EIGEN_SORT_ABS_ASC);

	for(int i = 0; i < 3; i++)
		eigenv[i] = gsl_vector_get(gsl_eigenvalues, i);

	complex<double> U_matter[3][3];
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			U_matter[i][j] = GSL_REAL(gsl_matrix_complex_get(gsl_U_matter, i, j)) + I*GSL_IMAG(gsl_matrix_complex_get(gsl_U_matter, i, j));
		}
	}

	/*freeeing GSL machinary*/
	gsl_eigen_hermv_free(workspace);
	gsl_matrix_complex_free(gsl_H);
	gsl_matrix_complex_free(gsl_U_matter);
	gsl_vector_free(gsl_eigenvalues);
	
	//printa_matriz(U);
	//printa_matriz(U_matter);

	nu_state nu;
	for(int i = 0; i < 3; i++)
		nu.index[i] = U_matter[0][i];

	return(nu);
}

double survival_prob_sun_general_adiabatic(double Rfrac, double E, osc_param op){
	/***************************************/
	/*EVALUATING THE AVERAGED SURVIVAL PROB*/
	/***************************************/
	double eingenv[3];

	/*Adiabatic evolution*/
	nu_state nu = sun_adiabatic_decoherent_evolution(Rfrac, E, op, eingenv);

	/*Incoherent average survival probability*/
	double resultado = 0;
	for(int i = 0; i < 3; i++)
		resultado += norm(op.U[0][i])*norm(nu.index[i]);

	return(resultado);
}

//double survival_prob_sun_adiabatic_Zplike(double Rfrac, double E, double V_e, double V_mu, double V_tau){
//	/***************************************/
//	/*EVALUATING THE AVERAGED SURVIVAL PROB*/
//	/***************************************/
//	complex<double> U[3][3];
//	seta_PMNS(U, false);
//	/*Matter potential parameters*/
//	double cte = 280.5147710743851;
//	double me = 0.511;
//	double mZprime = 1;
//	V_e *= cte / (2*E*me + mZprime*mZprime);
//	V_mu *= cte / (2*E*me + mZprime*mZprime);
//	V_tau *= cte / (2*E*me + mZprime*mZprime);
//	double V_NP[3][3] = {{V_e, 0, 0},
//						 {0, V_mu, 0},
//						 {0, 0, V_tau}};
//	double eingenv[3];
//	/*Adiabatic evolution*/
//	nu_state nu = sun_adiabatic_decoherent_evolution(Rfrac, E, V_NP, eingenv);
//	/*Incoherent average survival probability*/
//	double resultado = 0;
//	for(int i = 0; i < 3; i++)
//		resultado += norm(U[0][i])*norm(nu.index[i]);
//	return(resultado);
//}

//double masses_values(int index, double Rfrac, double E, double V_e, double V_mu, double V_tau){
//	/***************************************/
//	/*EVALUATING THE AVERAGED SURVIVAL PROB*/
//	/***************************************/
//	complex<double> U[3][3];
//	seta_PMNS(U, false);
//	/*Matter potential parameters*/
//	double V_NP[3][3] = {{V_e, 0, 0},
//						 {0, V_mu, 0},
//						 {0, 0, V_tau}};
//	double eigenv[3];
//	/*Adiabatic evolution*/
//	sun_adiabatic_decoherent_evolution(Rfrac, E, V_NP, eigenv);
//	return(eigenv[index]);
//}