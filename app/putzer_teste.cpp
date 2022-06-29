#include<iostream>
#include"putzer_algorithm.h"

int main(){
	complex<double> A[3][3] = {{-4.53146, -3.11439, -0.816852},
					  		   {-3.11439, -2.45472, -0.561407},
					  		   {-0.816852, -0.561407, -10.1369}};

	complex<double> w[3];

	complex<double> e1[3] = {1, 0, 0};
	complex<double> e2[3] = {0, 1, 0};
	complex<double> e3[3] = {0, 0, 1};

	
	matrixexp_times_vector_3x3_putzer(A, e1, w);
	for(int i = 0; i < 3; i++)
		cout << w[i] << " ";
	cout << endl;

	matrixexp_times_vector_3x3_putzer(A, e2, w);
	for(int i = 0; i < 3; i++)
		cout << w[i] << " ";
	cout << endl;

	matrixexp_times_vector_3x3_putzer(A, e3, w);
	for(int i = 0; i < 3; i++)
		cout << w[i] << " ";
	cout << endl;

	cout << endl;

	double lambda[3] = {1, 7, 3};

	if(lambda[0] > lambda[2])
		swap(lambda[0], lambda[2]);
	if(lambda[0] > lambda[1])
		swap(lambda[0], lambda[1]);
	if(lambda[1] > lambda[2])
		swap(lambda[1], lambda[2]);

	for(int i = 0; i < 3; i++)
		cout << lambda[i] << " ";
	cout << endl;
}