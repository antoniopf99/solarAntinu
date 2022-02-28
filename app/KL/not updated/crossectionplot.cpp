#include<cmath>
#include<fstream>
#include"crosssections.h"

using namespace std;

int main(){

	/***********************/
	/*OUTPUT FILE DIRECTORY*/
	/***********************/
	string outdir = "results/";

	/**************/
	/*NUMERIC PACE*/
	/**************/
	long double PASSO = 0.01;

	/****************************/
	/*MINIMUM FLUX ENERGRY (MeV)*/
	/****************************/
	 long double EMIN = 1.806; /*10^3.257 keV ~ 1.8 MeV inverse beta decay treshold*/

	/***************************/
	/*MAXIMUM FLUX ENERGY (MeV)*/
	/***************************/
	 long double EMAX = 200;

	ofstream saida;
	string nomesaida = "crosssection.txt";
 	saida.open(outdir + nomesaida);
 	for(long double E = EMIN; E < EMAX; E += PASSO)
 		saida << (E) << " " << inverseBeta(E) << endl;
 	saida.close();
 	saida.clear();
 	return(0);
}
