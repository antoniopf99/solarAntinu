#include<iostream>
#include<cmath>
#include<fstream>
#include"crosssections.h"
#include"reactors.h"
#include"oscilation.h"
#include"detectors.h"

void zBins(double Ntheo[][3]){

 	/*Number of bins*/
 	int numberOfBins = 17;

 	zeraBins(Ntheo, numberOfBins);
}

void calcBins(vector<Reactor> &reatores, double Ntheo[][3], bool oscila, long double theta12, long double theta13, long double Dmq){

	/*****************/
 	/*BINS PARAMETERS*/
 	/*****************/

 	/*Number of numeric intervals (MUST BE EVEN)*/
 	int numericIntervals = 500;

 	/*Number of bins*/
 	int numberOfBins = 17;
 	/*Size of the bins*/
 	double BinSize = 0.425;

 	/*Minimum energy considered in the analysis*/
 	double EanalMin =  1.6823; /*0.9 + 0.7823*/
 	/*Maximum energy considered in the analysis*/
 	double EanalMax = 15;

 	/*Name of the file with detector efficiency data*/
 	string datafile = "effkamland.txt";
 	/*Maximum distance from the reactor considered in the analysis, in cm*/
 	long double distance = 1000 * pow(10,5);
 	/*Resolution of the detector in %/sqrt(MeV)*/
 	double res = 6.5;

 	calculaBins(reatores, distance, res, datafile, numericIntervals, Ntheo, BinSize, numberOfBins, EanalMin, EanalMax, oscila, theta12, theta13, Dmq);
}

void adicionaBG(double Ntheo[][3], int numberOfBins){

	/*Name of the file with detector background data*/
 	string datafile = "data/BG.txt";

	double bg;
	ifstream arquivo;
	arquivo.open(datafile);
	for(int i = 0; i < numberOfBins; i++){
		arquivo >> bg;
		Ntheo[i][0] += bg;
	}
}

int main(){

	/***********************/	
	/*OUTPUT FILE DIRECTORY*/
	/***********************/
	string outdir = "results/";
	
	geoCoord ponto;
	/*KamLAND coordinates*/
	ponto.lat = 36.4225;
	ponto.lon = 137.3153;

	/*****************/
	/*BEST FIT VALUES*/
	/*****************/
	long double DmqBEST = 7.58 * pow(10,-5);
	long double theta12BEST = 0.642432;
	long double theta13BEST = 0.141897;

	int mais, maisanos, oscila, bg, regions;
	string nomesaida;
	ofstream saida;
	vector<Reactor> reatores;
	long double distancia;
	double Ntheo[17][3];
 	do{
 		zBins(Ntheo);
 		cout << "Enter with the output file name: ";
 		cin >> nomesaida;
 		saida.open(outdir + nomesaida);
		cout << "Do you want to consider oscilation effects? (1=y, 0=n): ";
 		cin >> oscila;
 		cout << endl;
 		cout << "Enter with the number of years that the experiment collected data: ";
 		cin >> maisanos;
 		for(int i = 0; i < maisanos; i++)
 			leOsDados(reatores, maisanos);
		for(int i = 0; (unsigned)i < reatores.size(); i++){
	 		distancia = Distancia(ponto,reatores[i].loc);
	 		reatores[i].distance = distancia;
 		}
 		cout << "Do you want a regions plot? (2 = tan2t12 vs Dmq, 1 = tan2t12 vs sin2t13, 0=n): ";
 		cin >> regions;
 		if(regions){
			/**********************************/
			/*   KamLAND DATA FOR EACH BIN    */
			/**********************************/
			double Nobs[17] = {143.17523, 230.09770, 229.07902, 238.79564, 170.27673, 196.16984, 189.26891, 177.36781, 182.23149, 123.12435, 77.25280, 70.35182, 49.33290, 15.52012, 10.38385, 3.18858, 2.31701};
			double d[17] = {11.76455, 15.00001, 15.00001, 15.73536, 12.64713, 13.67654, 13.52937, 12.94113, 13.67648, 11.02948, 8.82354, 8.23524, 7.05888, 3.82348, 3.82348, 3.82348, 3.82348};

			/*****************************/
			/*DETERMINATION OF THE REGION*/
			/*****************************/
			/*Delta m squared will be *e-4*/
			long double Dmqini = 0.2;
			long double Dmqfim = 1.2;
			/*theta12 scale is in tan^2*/
			long double theta12ini = 0.1;
			long double theta12fim = 1.1;
			/*theta13 scale is in sin^2*/
			long double theta13ini = 0;
			long double theta13fim = 0.8;
				
			/*RESOLUTION PACE OF THE REGION*/
			long double PASSO = 0.01;

			bool primeiro = true;
			if(regions == 2){
				float progresso = 0, progressotodo = ((theta12fim-theta12ini)/PASSO + 1)*((Dmqfim-Dmqini)/PASSO + 1);
				long double menor, tanthetamenor, Dmqmenor, tantheta, theta, cofDmq, Dmq, chi2;
				for(tantheta = theta12ini; tantheta <= theta12fim; tantheta += PASSO){
					if(tantheta != theta12ini)
						saida << endl;
					theta = atan(tantheta);
					for(cofDmq = Dmqini; cofDmq <= Dmqfim; cofDmq += PASSO){
						Dmq = cofDmq * pow(10,-4);
						zBins(Ntheo);
						calcBins(reatores, Ntheo, (bool)oscila, theta, theta13BEST, Dmq);
						adicionaBG(Ntheo, 17);
						chi2 = 0;
						for(int i = 0; i < 17; i++)
							chi2 += pow(Ntheo[i][0] - Nobs[i],2) / pow(d[i],2);
						if(primeiro || chi2 < menor){
							menor = chi2;
							tanthetamenor = tantheta;
							Dmqmenor = Dmq;
							if(primeiro)
								primeiro = false;
						}
						saida << pow(tantheta,2) << " " << Dmq << " " << chi2 << endl;
						progresso++;
						printf("%0.5f%% done.\n", progresso/progressotodo*100);
					}
				}
				cout << endl;
				cout << "Best fit " << "(chi^2 = " << menor << ") at tan^2theta = " << pow(tanthetamenor,2) << " and Delta m^2 = " << Dmqmenor << "." << endl;
			}
			if(regions == 1){
				float progresso = 0, progressotodo = ((theta12fim-theta12ini)/PASSO + 1)*((Dmqfim-Dmqini)/PASSO + 1);
				long double menor, tantheta12menor, sintheta13menor, tantheta12, theta12, sintheta13, theta13, chi2;
				for(tantheta12 = theta12ini; tantheta12 <= theta12fim; tantheta12 += PASSO){
					if(tantheta12 != theta12ini)
						saida << endl;
					theta12 = atan(tantheta12);
					for(sintheta13 = theta13ini; sintheta13 <= theta13fim; sintheta13 += PASSO){
						theta13 = asin(sintheta13);
						zBins(Ntheo);
						calcBins(reatores, Ntheo, (bool)oscila, theta12, theta13, DmqBEST);
						adicionaBG(Ntheo, 17);
						chi2 = 0;
						for(int i = 0; i < 17; i++)
							chi2 += pow(Ntheo[i][0] - Nobs[i],2) / pow(d[i],2);
						if(primeiro || chi2 < menor){
							menor = chi2;
							tantheta12menor = tantheta12;
							sintheta13menor = sintheta13;
							if(primeiro)
								primeiro = false;
						}
						saida << pow(tantheta12,2) << " " << pow(sintheta13,2) << " " << chi2 << endl;
						progresso++;
						printf("%0.5f%% done.\n", progresso/progressotodo*100);
					}
				}
				cout << endl;
				cout << "Best fit " << "(chi^2 = " << menor << ") at tan^2theta12 = " << pow(tantheta12menor,2) << " and sin^2theta13 = " << pow(sintheta13menor,2) << "." << endl;
			}
		}
		else{
			calcBins(reatores, Ntheo, (bool)oscila, theta12BEST, theta13BEST, DmqBEST);
 			reatores.clear();
 			cout << "Do you want to consider background effects? (1=y, 0=n): ";
 			cin >> bg;
 			if(bg)
 				adicionaBG(Ntheo, 17);
 			for(int i = 0; i < 17; i++){
 				saida << (Ntheo[i][1] - 0.7823) << " " << Ntheo[i][0] << endl;
 				saida << (Ntheo[i][2] - 0.7823) << " " << Ntheo[i][0] << endl;
 			}
 		}
 		saida.close();
 		saida.clear();
		cout << "Do you want to colect more data? (1=y, 0=n): ";
 		cin >> mais;
 		cout << endl;
 	}while(mais);
 	return(0);
}