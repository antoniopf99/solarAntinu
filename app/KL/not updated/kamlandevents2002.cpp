#include<iostream>
#include<cmath>
#include<fstream>
#include"crosssections.h"
#include"reactors.h"
#include"oscilation.h"
#include"detectors.h"

void zBins(double Ntheo[][3]){

 	/*Number of bins*/
 	int numberOfBins = 19;

 	zeraBins(Ntheo, numberOfBins);
}

void calcBins(vector<Reactor> &reatores, double Ntheo[][3], bool oscila, long double theta, long double Dmq, bool bestfit, double fraction){

	/*****************/
 	/*BINS PARAMETERS*/
 	/*****************/

 	/*Number of numeric intervals (MUST BE EVEN)*/
 	int numericIntervals = 1000;

 	/*Number of bins*/
 	int numberOfBins = 19;
 	/*Size of the bins*/
 	double BinSize = 0.425;

 	/*Minimum energy considered in the analysis*/
 	double EanalMin = 0.8323; /*0.05+0.7823*/
 	/*Maximum energy considered in the analysis*/
 	double EanalMax = 15;

 	/*Name of the file with detector efficiency data*/
 	string datafile = "";
 	/*Maximum distance from the reactor considered in the analysis, in cm*/
 	long double distance = 35000000;
 	/*Resolution of the detector in %/sqrt(MeV)*/
 	double res = 7.5;

 	calculaBins(reatores, distance, fraction, res, datafile, numericIntervals, Ntheo, BinSize, numberOfBins, EanalMin, EanalMax, oscila, theta, Dmq, bestfit);
}

int main(){

	/***********************/	
	/*OUTPUT FILE DIRECTORY*/
	/***********************/
	string outdir = "results/";
	
	/*********************************/
	/*PRECISION (# OF SUBINTERVALS   */
	/*IN EACH EMIN - EMAX INTERVAL)  */
	/*********************************/
	long N = 500;/*must be even!*/
	N *= 2;

	int mais, oscila, bins, regions;
	string nomesaida;
	ofstream saida;
	vector<Reactor> reatores;
	long double events, distancia;
	leOsDados(reatores);
 	cout << "File read successfully." << endl;
 	geoCoord ponto;

	/*KamLAND coordinates*/
	ponto.lat = 36.4225;
	ponto.lon = 137.3153;


	for(int i = 0; (unsigned)i < reatores.size(); i++){
	 	distancia = Distancia(ponto,reatores[i].loc);
	 	reatores[i].distance = distancia;
 	}
 	do{
 		cout << "Enter with the output file name: ";
 		cin >> nomesaida;
 		saida.open(outdir + nomesaida);
		cout << "Continuous (0) or bined (1) spectrum? " << endl;
 		cin >> bins;
		cout << "Do you want to consider oscilation effects? (1=y, 0=n): ";
 		cin >> oscila;
 		if(oscila && bins){
 			cout << "Do you want to construct an allowed regions map? (1=y, 0=n): ";
 			cin >> regions;
 		}
 		else
 			regions = 0;
 		cout << endl;
 		if(bins){
			double Ntheo[19][3];
			if(regions){
				/**********************************/
				/*   KamLAND DATA FOR EACH BIN    */
				/*(data at KamLAND start at bin 7)*/
				/**********************************/
				double Nobs[19] = {0,0,0,0,0,0,6.5,11,8.5,7.5,7.5,4,4.6,2.1,0,0,0,0,0};

				/*****************************/
				/*DETERMINATION OF THE REGION*/
				/*****************************/
				/*Delta m squared will be in log scale*/
				long double Dmqini = -6;
				long double Dmqfim = 0;
				/*theta scale is in tan^2*/
				long double thetaini = 0;
				long double thetafim = 1.1;
				
				/*RESOLUTION PACE OF THE REGION*/
				long double PASSO = 0.01;

				float progresso = 0, progressotodo = ((thetafim-thetaini)/PASSO + 1)*((Dmqfim-Dmqini)/PASSO + 1);
				long double menor, tanthetamenor, Dmqmenor, tantheta, theta, pot, Dmq, chi2, chiG2, chiP2, d2[19];
				bool primeiro = true;
				for(tantheta = thetaini; tantheta <= thetafim; tantheta += PASSO){
					if(tantheta != thetaini)
						saida << endl;
					theta = atan(tantheta);
					for(pot = Dmqini; pot <= Dmqfim; pot += PASSO){
						Dmq = pow(10,pot);
						zBins(Ntheo);
						calcBins(reatores, Ntheo, (bool)oscila, theta, Dmq, false, 0.3975);
						for(int i = 6; i < 19; i++)
							d2[i] = Nobs[i] + pow(Nobs[i]*0.0642,2);
						chiG2 = 0;
						chiP2 = 0;
						for(int i = 6; i < 19; i++){
							if(4 <= Nobs[i])
								chiG2 += pow(Nobs[i] - Ntheo[i][0],2) / d2[i];
							else{
								if(Nobs[i] > 0)
									chiP2 += 2*(Ntheo[i][0] - Nobs[i]) + 2*Nobs[i]*log(Nobs[i]/Ntheo[i][0]);
								else
									chiP2 += 2*(Ntheo[i][0] - Nobs[i]);
							}
						}
						chi2 = chiG2 + chiP2;
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
				cout << "Best fit " << "(chi^2 = " << menor << ") at tan^2theta = " << tanthetamenor << " and Delta m^2 = " << Dmqmenor << "." << endl;
			}
			else{
				zBins(Ntheo);
 				calcBins(reatores, Ntheo, (bool)oscila, 0, 0, true, 0.3975);
 				for(int i = 0; i < 19; i++){
 					saida << (Ntheo[i][1] - 0.7823) << " " << Ntheo[i][0] << endl;
 					saida << (Ntheo[i][2] - 0.7823) << " " << Ntheo[i][0] << endl;
 				}
 			}
 		}
 		else{
 			/****************************/
			/*MINIMUM FLUX ENERGRY (MeV)*/
			/****************************/
			double EMIN = 0; /*10^3.257 keV ~ 1.8 MeV inverse beta decay treshold*/
			
			/***************************/
			/*MAXIMUM FLUX ENERGY (MeV)*/
			/***************************/
			double EMAX = 10;

 			long double h = (EMAX - EMIN) / N;
 			for(long double E = EMIN; E < EMAX; E += h){
 				events = 0;
 				for(Reactor it : reatores){
 					if(!it.country.compare("JP")){
 						if(oscila)
 							events += exposureXsection(E, 0.3975) * survivalProb1GEN(it.distance * pow(10,-5)/*cm to km*/,E * pow(10,-3)/*MeV to GeV*/,0,0,true) * it.Fluxo(E);
 						else
 							events += exposureXsection(E, 0.3975) * it.Fluxo(E);
 					}
 				}
 				saida << (E - 0.7823) << " " << events << endl;
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