#ifndef PUTZERT_ALGORITHM_H
#define PUTZERT_ALGORITHM_H

#include<iostream>
#include<cmath>
#include<complex>

using namespace std;

#ifndef IMAGINARY_UNIT
#define IMAGINARY_UNIT
const complex<double> I = complex<double>(0.0, 1.0);
#endif //IMAGINARY_UNIT


/*
Computes exponential of the matrix A times v, and returns the result in result.
result = e^A v
*/
void matrixexp_times_vector_3x3_putzer(complex<double> A[][3], complex<double> v[3], complex<double> result[3]);

#endif