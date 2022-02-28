#include<iostream>
#include<cmath>
#include<fstream>
#include"crosssections.h"
#include"reactors.h"
#include"oscilation.h"
#include"detectors.h"

using namespace std;

long double taishan1Flux(long double energy){
	int type = 0; /*PWR reactor*/
	int MOX = 0; /*reactor doesn't use MOX*/
	long double thermalPower = 4590; /*in MW*/
	long double loadFactor = 100; /*assuming the reactor will be fully loaded at all times*/
	long double distance = 3000; /*3000cm, 30m*/
	return(customFlux(type, MOX, thermalPower, loadFactor, energy, distance));
}

long double taishan2Flux(long double energy){
	int type = 0; /*PWR reactor*/
	int MOX = 0; /*reactor doesn't use MOX*/
	long double thermalPower = 4590; /*in MW*/
	long double loadFactor = 100; /*assuming the reactor will be fully loaded at all times*/
	long double distance = 16000; /*16000cm, 160m*/
	return(customFlux(type, MOX, thermalPower, loadFactor, energy, distance));
}

int main(){

	/***********************/	
	/*OUTPUT FILE DIRECTORY*/
	/***********************/
	string outdir = "results/";

	/*energy range in MeV where the flux will be calculated*/
	double minE = 0;
	double maxE = 10;

	/*numeric step*/
	double h = 0.001;

	ofstream saida;

	double spec;

	///*only taishan-1*/
	//string t1 = "t1.txt";
	//saida.open(outdir + t1);
	//for(double E = minE; E < maxE; E += h){
	//	flux = taishan1Flux(E); 
	//	saida << E << " " << flux << endl;
	//}
	//saida.close();
	///*taishan-1 and 2*/
	//string t2 = "t1andt2.txt";
	//saida.open(outdir + t2);
	//for(double E = minE; E < maxE; E += h){
	//	flux = taishan1Flux(E) + taishan2Flux(E); 
	//	saida << E << " " << flux << endl;
	//}
	//saida.close();

	/*taishan-1 and 2*/
	string spectrum = "taospectrum.txt";
	saida.open(outdir + spectrum);
	/*TAO has 7.2e28 protons, the cross section is missing a factor of 1e-42*/
	long double Np = 7.2 * pow(10,-14);
	for(double E = minE; E < maxE; E += h){
		spec = Np * taishan1Flux(E) * inverseBeta(E); 
		saida << E << " " << spec << endl;
	}
	saida.close();

	return(0);
}