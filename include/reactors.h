#ifndef REACTORS_H
#define REACTORS_H

#include<iostream>
#include<fstream>
#include<vector>
#include<cstring>
#include<string>
#include"geoCoords.h"
#include"energyConverts.h"

using namespace std;

/*REACTOR CLASS*/
class Reactor{
	public:
		geoCoord loc;

		/*
		Load factor in [0,1] (averaged over the year),
		thermal power in MW,
		distance in cm. 
		*/
		long double loadFactor, thermalPower, distance;
		int type, MOX;
		string country, name, tipo;

		/*RETURNS DIFFERENTIAL FLUX dF/dE IN 1/cm² MeV s, ENERGY IN MeV*/
		long double Fluxo(long double E);
};

/*LAMBDA FUNCTION: (ISOTOPE, ENERGY IN keV)*/
long double Lambda(int i, long double E);

/*FUNCTION TO MAKE THE FLUX CALCULATION FUNCTION MORE CONCISE: (REACTOR TYPE, ENERGY IN keV)*/
long double FatorA(int i, long double E);

/*quantosAnos IS HOW MANY YEARS THE DETECTOR RAN, TO CALCULATE THE FRACTIONAL CONTRIBUTION OF A SPECIFIC YEAR,
SINCE EACH DATABASE READ CORRESPONDS TO A SPECIFIC YEAR*/
void leOsDados(vector<Reactor> &reatores, int quantosAnos);

/*CALCULATES THE FLUX FROM A SINGLE CUSTOM REACTOR*/
long double customFlux(int type, int MOX, long double thermalPower, long double loadFactor, long double E, long double distance);

#endif
