#ifndef CROSSSECTIONS_H
#define CROSSSECTIONS_H

#include<cmath>

//typedef struct XSECTION_PARAM{
//	double e;
//	double mZprime;
//	double 
//} xsection_param;

typedef struct ZPRIME_PARAM{
	bool antinu;
	bool NC_and_CC;
	double Enu;
	double T;
	double gX;
	double a;
	double e;
	double mZprime;
	double QlL;
	double QlR;
	double QnuL;
} zprime_param;

/*Returns the cross section of the inverse beta decay times 10^42,
as a function of the incoming neutrino energy in MeV, in units of cm^2*/
double inverseBeta(double E);

double nueScatteringJUNO(zprime_param p);

/*Returns the cross section of antinue and e scattering,
as a function of the incoming neutrino energy Enu in MeV and the
electron recoil energy T in MeV, in units of MeV^-1cm^2*/
double nueScatteringSM(double Enu, double T);

/*Returns the cross section of antinue and e scattering times 10^45,
as a function of the incoming neutrino energy Enu in MeV and the
electron recoil energy T in MeV, in units of MeV^-1cm^2*/
double nueScatteringBSM(zprime_param p);

#endif