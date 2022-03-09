#ifndef MAGNUS_SOLAR_H
#define MAGNUS_SOLAR_H

#include<cmath>
#include<complex>
#include"solarmodel.h"

using namespace std;

const complex<double> I = complex<double>(0.0,1.0);

double survival_prob_sun_magnus(double Rfrac, double E, double V_e, double V_mu, double V_tau);

#endif
