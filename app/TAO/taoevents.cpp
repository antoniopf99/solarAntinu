#include<iostream>
#include<cmath>
#include<fstream>
#include"crosssections.h"
#include"reactors.h"
#include"oscilation.h"
#include"detectors.h"

using namespace std;

/*******************/
/*NUMERIC PRECISION*/
/*******************/
int N = 100;


//long double numInt(long double a, long double b, long double (*f)(long double, long double, long double)){
//	long double soma, h;	
//	h = (b-a) / N;
//	soma = f(a) + f(a + N*h);
//	for(int i = 1; i < N; i+=2)
//		soma += 4*f(a + i*h) + 2*f(a + (i+1)*h);
//	soma *= h/3;
//	return(soma);
//}

long double taishan1Flux(long double energy){
	int type = 0; /*PWR reactor*/
	int MOX = 0; /*reactor doesn't use MOX*/
	long double thermalPower = 4590; /*in MW*/
	long double loadFactor = 100; /*assuming the reactor will be fully loaded at all times*/
	long double distance = 3000; /*3000cm, 30m*/
	return(customFlux(type, MOX, thermalPower, loadFactor, energy, distance));
}

long double partialX(long double Enu, long double T1, long double T2){
	/*TAO has 7.2e28 protons, a day is 8.64e4s, to together they form 6.2208e23,
	the cross section is missing a factor of 1e-45. 80% eff*/
	long double Np = 0.8 * 6.2208 * pow(10,-14);
	
	long double soma, h;	
	h = (T2-T1) / N;
	soma = antinueScatteringSM(Enu, T1) + antinueScatteringSM(Enu, T1 + N*h);
	for(int i = 1; i < N; i+=2)
		soma += 4*antinueScatteringSM(Enu, T1 + i*h) + 2*antinueScatteringSM(Enu, T1 + (i+1)*h);
	soma *= (h/3)*Np;
	return(soma);
}

int main(){

	/***********************/	
	/*OUTPUT FILE DIRECTORY*/
	/***********************/
	string outdir = "results/";

	/*energy range in MeV where the flux will be calculated*/
	double minE = 0;
	double maxE = 10;

	/*energies involved in the analysis*/
	double EAnalMin = 0;
	double EAnalMax = 20;

	/*number of bins*/
	int nBins = 20;

	double binSize = (maxE - minE) / nBins;

	ofstream saida;

	string spectrum = "taoevents.txt";
	saida.open(outdir + spectrum);

	long double soma, h, T1, T2;
	h =  (EAnalMax - EAnalMin)/ N;
	for(int i = 0; i < nBins; i++){
		T1 = minE + i*binSize;
		T2 = minE + (i + 1)*binSize;
		soma = taishan1Flux(EAnalMin)*partialX(EAnalMin, T1, T2) + taishan1Flux(EAnalMin + N*h)*partialX(EAnalMin + N*h, T1, T2);
		for(int j = 1; j < N; j+=2)
			soma += 4*taishan1Flux(EAnalMin + j*h)*partialX(EAnalMin + j*h, T1, T2) + 2*taishan1Flux(EAnalMin + (j+1)*h)*partialX(EAnalMin + (j+1)*h, T1, T2);
		soma *= h/3;
		saida << T1 << " " << soma << endl;
		saida << T2 << " " << soma << endl;
	}
	saida.close();

	return(0);
}