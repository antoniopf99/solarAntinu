#include<iostream>
#include<cmath>
#include<fstream>
#include"oscilation.h"

using namespace std;

int main(){

	/***********************/	
	/*OUTPUT FILE DIRECTORY*/
	/***********************/
	//string file, outdir = "results/";
	//ofstream saida;

	seta_dados_Ne();

	double L, E, *p = NULL;

	L = 696342;

	E = 0.01;

	string file = "sunnumeric10MeV.txt";

	cout << "A prob de sobrevivência deu " << survival_prob_sun(E, 0, L, 0.1, p, true, file) << endl;

	//file = "sunnumeric.txt";
	//saida.open(outdir + file);
	//for(double E = 0.0001; E < 1; E *= 10)
	//	cout << "Resultado para " << E*1000 << " MeV: " << survival_prob_sun(E, 0, L, 0.1, p, false, "") << endl;
	//double E;
	//file = "oscilaneutrino1GeV.txt";
	//saida.open(outdir + file);
	//for(double pot = -1; pot < 1.5; pot += 0.01){
	//	E = pow(10,pot);
	//	saida << E << " " << survivalProb_3genSun_Ad(E, 0) << endl; 
	//}
	//saida.close();
	//file = "boronaverage48.txt";
	//saida.open(outdir + file);
	//	for(double pot = -1; pot < 1.5; pot += 0.01){
	//	E = pow(10,pot);
	//	saida << E << " " << survivalProb_3genSun_Ad_Boron(E) << endl; 
	//}
	//saida.close();
	
	return(0);
}