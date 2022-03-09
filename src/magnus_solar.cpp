#include<iostream>
#include"magnus_solar.h"

double survival_prob_sun_magnus(double Rfrac, double E, double V_e, double V_mu, double V_tau){
	/************************/
	/*PMNS MATRIX PARAMETERS*/
	/************************/
	double theta12 = 33.45 * M_PI / 180;
	double theta13 = 8.62 * M_PI / 180;
	double theta23 = 42.1 * M_PI / 180;
	double deltaCP = 195.0 * M_PI / 180;

	double s12 = sin(theta12);
	double s13 = sin(theta13);
	double s23 = sin(theta23);
	double c12 = cos(theta12);
	double c13 = cos(theta13);
	double c23 = cos(theta23);

	/************************************/
	/*MASS-SQUARED DIFFERENCES (IN eV^2)*/
	/************************************/
	double Dmq21 = 7.42e-5;
	double Dmq31 = 2.514e-3;

	/*******************************/
	/*PMNS MATRIX, SCCORDING TO PDG*/
	/*******************************/
	complex<double> U[3][3], faseCP = exp(I*deltaCP);

	U[0][0] = c12 * c13;
	U[0][1] = s12 * c13;
	U[0][2] = s13 * conj(faseCP);
	U[1][0] = -s12 * c23 - c12 * s13 * s23 * faseCP;
	U[1][1] = c12 * c23 - s12 * s13 * s23 * faseCP;
	U[1][2] = c13 * s23;
	U[2][0] = s12 * s23 - c12 * s13 * c23 * faseCP;
	U[2][1] = -c12 * s23 - s12 * s13 * c23 * faseCP;
	U[2][2] = c13 * c23;
	
	/****************************/
	/*SETTING UP THE HAMILTONIAN*/
	/*(in units of solar radius)*/
	/****************************/
	
	/*the vacuum portion H_0*/
	double H_0[3][3];
	double cte = 1.7644318e9;
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			H_0[i][j] = 0.0;
	H_0[1][1] = cte*Dmq21/E;
	H_0[2][2] = cte*Dmq31/E;

	/*the matter potential, in units of sqrt(2)*G_F*n_e*/
	complex<double> W[3][3];
	double V[3][3];
	double cte2 = 2.6933932e2;
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			W[i][j] = V[i][j] = 0.0;

	V[0][0] = cte2*V_e;
	V[1][1] = cte2*V_mu;
	V[2][2] = cte2*V_tau;

	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			for(int k = 0; k < 3; k++) 
					W[i][j] += conj(U[k][i])*V[k][k]*U[k][j]; /*assuming that V is diagonal!*/

	//for(int i = 0; i < 3; i++){
	//	for(int j = 0; j < 3; j++){
	//		cout << W[i][j] - conj(W[j][i]) << " ";
	//	}
	//	cout << endl;
	//}

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
		Vatual[i] = U[0][i]; /*0, electron neutrino!*/

	/****************************/
	/*PERFORMING THE INTEGRATION*/
	/****************************/
	complex<double> Omega[3][3];
	double xi_plus, xi_minus, v_plus, v_minus, h = 0.00005, sqrt3 = sqrt(3);

	for(double xi = Rfrac; xi < 1; xi += h){
		for(int i = 0; i < 3; i++)
			for(int j = 0; j < 3; j++)
				Omega[i][j] = 0.0;
		
		/****************************/
		/*COMPUTING THE OMEGA MATRIX*/
		/****************************/
		xi_plus = xi + (1.0 + (1.0/sqrt3))*(h/2.0);
		xi_minus = xi + (1.0 - (1.0/sqrt3))*(h/2.0);
		v_plus = pow(10.0, electron_density_sun_NA(xi_plus));
		v_minus = pow(10.0, electron_density_sun_NA(xi_minus));

		/*i include an extra -i factor to make it hermitian*/
		for(int i = 0; i < 3; i++){
			for(int j = 0; j < 3; j++){
				Omega[i][j] = -I*(-I*(H_0[i][j] + 0.5*(v_plus+v_minus)*W[i][j])*h + (sqrt3/12.0)*(v_plus-v_minus)*H_0cW[i][j]*h*h);
			}
		}

		//for(int i = 0; i < 3; i++){
		//	for(int j = 0; j < 3; j++){
		//		cout << Omega[i][j] - conj(Omega[j][i]) << " ";
		//	}
		//	cout << endl;
		//}
		//cout << endl;
	
		/**************************************/
		/*CALCULATING THE EXPONENTIAL OF OMEGA*/
		/**************************************/
		/*getting the trace*/
		complex<double> A_0[3][3], trace_Omega = 0;
		for(int i = 0; i < 3; i++){
			trace_Omega += Omega[i][i];
			for(int j = 0; j < 3; j++)
				A_0[i][j] = 0.0;
		}

		/*tracelss matrix A_0*/
		for (int i = 0; i < 3; i++){
			for(int j = 0; j < 3; j++){
				if(i == j)
					A_0[i][i] = Omega[i][i] - trace_Omega/3.0;
				else
					A_0[i][j] = Omega[i][j];
			}
		}

		/*parameters of the Putzer algorithm*/
		double lambda[3];
		double p = 0.0, q, sqrtpo3, arg_acos, acos_value; /*notice that those are real anyway!*/
		
		for(int i = 0; i < 3; i++)
			for(int j = 0; j < 3; j++)
				p += real(A_0[i][j]*A_0[j][i])/2;

		q = -real(A_0[0][0]*(A_0[1][1]*A_0[2][2] - A_0[1][2]*A_0[2][1]) - A_0[0][1]*(A_0[1][0]*A_0[2][2] - A_0[1][2]*A_0[2][0]) + A_0[0][2]*(A_0[1][0]*A_0[2][1] - A_0[1][1]*A_0[2][0])); 
		sqrtpo3 = sqrt(p/3.0);
		arg_acos = ((3.0*q)/(2.0*p))*(1.0/sqrtpo3);
		for(int i = 0; i < 3; i++){
			if(arg_acos < -1.0)
				acos_value = M_PI;
			else if(arg_acos > 1.0)
				acos_value = 0.0;
			else
				acos_value = acos(arg_acos);
			if(q <= 0.0)
				lambda[i] = 2.0*sqrtpo3*cos((1.0/3.0)*acos_value - 2.0*M_PI*i/3.0);
			else
				lambda[i] = -2.0*sqrtpo3*cos((1.0/3.0)*acos_value - 2.0*M_PI*i/3.0);
		}



		if(lambda[0] > lambda[2])
			swap(lambda[0], lambda[2]);
		if(lambda[0] > lambda[1])
			swap(lambda[0], lambda[1]);
		if(lambda[1] > lambda[2])
			swap(lambda[1], lambda[2]);

		double a = lambda[1] - lambda[0];
		double b = lambda[2] - lambda[0];
		complex<double> r0 = -(1.0 - exp(I*a))/a;
		complex<double> r1 = -(1.0/(a - b))*((1.0 - exp(I*a))/a - (1.0 - exp(I*b))/b);

		/*computing the next step*/
		complex<double> expi_l0_z = exp(I*(lambda[0] + trace_Omega/3.0));
		complex<double> A_0Vatual[3], exp_omega[3][3], Id_mat[3][3];
		for(int i = 0; i < 3; i++){
			Vprox[i] = A_0Vatual[i] = 0.0;
			for(int j = 0; j < 3; j++){
				exp_omega[i][j] = 0;
				if(i == j)
					Id_mat[i][j] = 1;
				else
					Id_mat[i][j] = 0;
			}
		}

		/*the matrix exponential*/
		//for(int i = 0; i < 3; i++){
		//	for(int j = 0; j < 3; j++){
		//		exp_omega[i][j] += exp(lambda)
		//		for(int k = 0; k < 3; k++)
		//					
		//	}
		//}

		///*********************************/
		//complex<double> teste_unitario[3][3];
		//for(int i = 0; i < 3; i++)
		//	for(int j = 0; j < 3; j++)
		//		teste_unitario[i][j] = 0;
		//for(int i = 0; i < 3; i++)
		//	for(int j = 0; j < 3; j++)
		//		for(int k = 0; k < 3; k++)
		//			teste_unitario[i][j] += exp_omega[i][k]*conj(exp_omega[j][k]);
		//for(int i = 0; i < 3; i++){
		//	for(int j = 0; j < 3; j++)
		//		cout << teste_unitario[i][j] << " ";
		//	cout << endl;
		//}
		//cout << endl;
		///*********************************/

		/*the next step*/
		for(int i = 0; i < 3; i++)
			for(int j = 0; j < 3; j++)
				Vprox[i] += exp_omega[i][j]*Vatual[j];

		for(int i = 0; i < 3; i++){
			Vprox[i] += expi_l0_z*(1.0 - lambda[0]*(r0 - lambda[1]*r1))*Vatual[i]; /*portion proportional to the identity*/
			for(int j = 0; j < 3; j++)
				A_0Vatual[i] += A_0[i][j]*Vatual[j];
		}
		for(int i = 0; i < 3; i++){ /*portion proportional to A_0*/
			Vprox[i] += expi_l0_z*(r0 + lambda[2]*r1)*A_0Vatual[i];
			for(int j = 0; j < 3; j++) /*portion proportional to A_0^2*/
				Vprox[i] += expi_l0_z*r1*A_0[i][j]*A_0Vatual[j];
		}
	
		/*setting back current step*/
		for(int i = 0; i < 3; i++)
			Vatual[i] = Vprox[i];
	}

	/***************************************/
	/*EVALUATING THE AVERAGED SURVIVAL PROB*/
	/***************************************/
	double resultado = 0;

	for(int i = 0; i < 3; i++)
		resultado += norm(U[0][i])*norm(Vatual[i]);

	return(resultado);
}