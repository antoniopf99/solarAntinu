#include<iostream>
#include<cmath>
#include<fstream>
#include"solarflux.h"
#include"crosssections.h"
#include"detectors.h"
#include"oscilation.h"

using namespace std;

long double Xsection(long double Enu, long double T, double *p){
	//return(1);
	return(antinueScatteringSM(Enu, T));
}

long double survivalProb(long double Enu, double *p){
	return(survivalProb_2genSun_Ad(Enu,0));
	//return(1);
}

long double Flux(long double Enu){
	return(solarFlux(3, Enu));
}

long double countRate(long double T){
	double *xp, *spp;
	return(countRate_nu_e_Scattering(&Flux, &Xsection, xp, &survivalProb, spp, T));
}

long double Ntau_func(long double E){
	double Nelectrons1 = 7.9*3.38*pow(10,32);
	double Nelectrons2 = 12.2*3.38*pow(10,32);
	double Nelectrons3 = 16.5*3.38*pow(10,32);
	double year = 3.154*pow(10,7);
	year *= 10;
	double Ntau1 = Nelectrons1*year;
	double Ntau2 = Nelectrons2*year;
	double Ntau3 = Nelectrons3*year;
	double Ntau;
	if(E <= 3)
		Ntau = Ntau1;
	else if(E <= 5)
		Ntau = Ntau2;
	else
		Ntau = Ntau3;
	return(Ntau);
}

int main(){
	/***********************/	
	/*OUTPUT FILE DIRECTORY*/
	/***********************/
	string outdir = "results/";

	ofstream saida;
	string file = "junoeventsosc.txt";
	saida.open(outdir + file);

	int numberOfBins = 140;
	double binMin = 2;
	double binMax = 16;
	double Bins[210];
	double Eres = 3;

	double binsize = (binMax-binMin)/numberOfBins;

	double *xp = NULL, *spp = NULL;

	calcula_eventos_Bins(Bins, binMin, binMax, numberOfBins, &Ntau_func, &Flux, &Xsection, xp, &survivalProb, spp, Eres);

	for(int i = 0; i < numberOfBins; i++){
		saida << binMin + i*binsize << " " << Bins[i] << endl;
		saida << binMin + (i+1)*binsize << " " << Bins[i] << endl;;
	}

	saida.close();
	return(0);
}