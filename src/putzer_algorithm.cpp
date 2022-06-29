#include"putzer_algorithm.h"

void matrixexp_times_vector_3x3_putzer(complex<double> A[][3], complex<double> v[3], complex<double> result[3]){
		/*constructing a traceless matrix A_0*/
		complex<double> trace_A = 0;
		complex<double> A_0[3][3];

		for(int i = 0; i < 3; i++){
			trace_A += A[i][i];
			for(int j = 0; j < 3; j++)
				A_0[i][j] = 0.0;
		}

		for (int i = 0; i < 3; i++){
			for(int j = 0; j < 3; j++){
				if(i == j)
					A_0[i][i] = A[i][i] - trace_A/3.0;
				else
					A_0[i][j] = A[i][j];
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
			lambda[i] = -2.0*sqrtpo3*cos((1.0/3.0)*acos_value - 2.0*M_PI*i/3.0);
			//if(q > 0)
			//	lambda[i] *= -1;
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
		complex<double> expi_l0_z = exp(I*(lambda[0] + trace_A/3.0));
		complex<double> A_0v[3], Id_mat[3][3];
		for(int i = 0; i < 3; i++){
			result[i] = A_0v[i] = 0.0;
			for(int j = 0; j < 3; j++){
				if(i == j)
					Id_mat[i][j] = 1;
				else
					Id_mat[i][j] = 0;
			}
		}

		/*computing the product*/
		for(int i = 0; i < 3; i++){
			result[i] += expi_l0_z*(1.0 - lambda[0]*(r0 - lambda[1]*r1))*v[i]; /*portion proportional to the identity*/
			for(int j = 0; j < 3; j++)
				A_0v[i] += A_0[i][j]*v[j];
		}
		for(int i = 0; i < 3; i++){ /*portion proportional to A_0*/
			result[i] += expi_l0_z*(r0 + lambda[2]*r1)*A_0v[i];
			for(int j = 0; j < 3; j++) /*portion proportional to A_0^2*/
				result[i] += expi_l0_z*r1*A_0[i][j]*A_0v[j];
		}
}