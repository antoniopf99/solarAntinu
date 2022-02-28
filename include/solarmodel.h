#ifndef SOLARMODEL_H
#define SOLARMODEL_H


#include<fstream>
#include<cmath>
#include<cstring>
#include<vector>

using namespace std;

typedef struct dat{
	double Rfrac;
	double dat;
	}data_point;

//	 reads file with data of solar electron number density
//	 DON'T CALL IT MULTIPLE TIMES
void seta_dados_Ne();

//   Sun electron number density in units of cm^{-3}/N_A, as a function 
//	 of the distance from the center of the sun in fracttion of sun's radius.
//	 Vector with data need to be provided, an interpolation will be performed.
//	 seta_dados_Ne() MUST BE CALLED
long double electron_density_sun_NA(long double Rfrac);

//   Logarithm (to the base 10) of the electron density in units of
//   cm^{-3}/N_A, where N_A is the Avogadro number
long double interpola_electron_density(long double Rfrac);

//ONLY WORKS UNTIL 0.5 OF SUNS RADIUS, REPRESENTET BY Rfrac
//0) Temperature in units of 10^6 deg (K)
//1) Logarithm (to the base 10) of the electron density in units of
//   cm^{-3}/N_A, where N_A is the Avogadro number
//2) Mass of the zone with the given radius in units of one solar mass 
//3) X(^7Be): beryllium-7 mass fraction
//   [For special purposes, it is useful to know the 7^Be mass fraction.]
//4) Fraction of pp neutrinos produced in the zone 
//5) Fraction of boron 8 neutrinos produced in the zone
//6) Fraction of nitrogen 13 neutrinos produced in the zone
//7) Fraction of oxygen 15 neutrinos produced in the zone
//8) Fraction of florine 17 neutrinos produced in the zone
//9) Fraction of beryllium 7 neutrinos produced in the zone
//10) Fraction of pep neutrinos produced in the zone
//11) Fraction of hep neutrinos produced in the zone
long double interpola_dados_solares(long double Rfrac, int dado);

#endif