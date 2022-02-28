#include<iostream>
#include<cmath>
#include<cstring>
#include<fstream>
#include"solarflux.h"

using namespace std;

int main(){

	/***********************/	
	/*OUTPUT FILE DIRECTORY*/
	/***********************/
	string outdir = "results/";

	double expmin = -1;
	double expmax = 2;
	double step = 0.01;
	double Enu;

	ofstream saida;

	string file = "fluxpp.txt";

	saida.open(outdir + file);

	for(double exp = expmin; exp <= expmax; exp += step){
		Enu = pow(10,exp);
		saida << Enu << " " <<  solarFlux(1, Enu) << endl;
	}

	saida.close();

	file = "fluxhep.txt";

	saida.open(outdir + file);

	for(double exp = expmin; exp <= expmax; exp += step){
		Enu = pow(10,exp);
		saida << Enu << " " << solarFlux(2, Enu) << endl;
	}

	saida.close();

	file = "flux8B.txt";

	saida.open(outdir + file);

	for(double exp = expmin; exp <= expmax; exp += step){
		Enu = pow(10,exp);
		saida << Enu << " " << solarFlux(3, Enu) << endl;
	}

	saida.close();

	file = "flux13N.txt";

	saida.open(outdir + file);

	for(double exp = expmin; exp <= expmax; exp += step){
		Enu = pow(10,exp);
		saida << Enu << " " << solarFlux(4, Enu) << endl;
	}

	saida.close();

	file = "flux15O.txt";

	saida.open(outdir + file);

	for(double exp = expmin; exp <= expmax; exp += step){
		Enu = pow(10,exp);
		saida << Enu << " " << solarFlux(5, Enu) << endl;
	}

	saida.close();

	file = "flux17F.txt";

	saida.open(outdir + file);

	for(double exp = expmin; exp <= expmax; exp += step){
		Enu = pow(10,exp);
		saida << Enu << " " << solarFlux(6, Enu) << endl;
	}

	saida.close();

	file = "fluxtotal.txt";

	saida.open(outdir + file);

	for(double exp = expmin; exp <= expmax; exp += step){
		Enu = pow(10,exp);
		saida << Enu << " " << solarFlux(0, Enu) << endl;
	}

	saida.close();

	return(0);
}