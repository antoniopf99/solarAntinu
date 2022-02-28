#include<iostream>
#include<cmath>
#include<fstream>
#include"reactors.h"
#include"oscilation.h"

int main(){

	/***********************/
	/*OUTPUT FILE DIRECTORY*/
	/***********************/
	string outdir = "results/";

	/**************/
	/*NUMERIC PACE*/
	/**************/
	long double PASSO = 0.001;

	/*************************************/
	/*MINIMUM FLUX ENERGRY, 10^PMIN (MeV)*/
	/*************************************/
	 long double PMIN = -1; /*10^2 keV Yuber*/

	/************************************/
	/*MAXIMUM FLUX ENERGY, 10^PMAX (KeV)*/
	/************************************/
	 long double PMAX = 1.6;/*10^4.6 keV Yuber*/

	/**************/
	/*MINIMUM FLUX*/
	/**************/
	 long double FMIN = pow(10,-20);

	string nomesaida;
	ofstream saida;
	vector<Reactor> reatores;
	long double x, flux, fluxparticular;
	leOsDados(reatores);
 	cout << "File read successfully." << endl;
	int mais = 1, oscila;
	do{
 		cout << "Enter with the output file name: ";
 		cin >> nomesaida;
 		saida.open(outdir + nomesaida);
 		flux  = 1;
 		long double distancia;
 		long double distanciamaisperto;
 		bool primeiro = true;
 		string nomemaisperto;
		geoCoord ponto;
 		cout << "Enter with the coordinates of the point where the flux will be calculated: ";
 		cin >> ponto.lat >> ponto.lon;
 		cout << "Do you want to consider oscilation effects? (1=y, 0=n): ";
 		cin >> oscila;
 		for(int i = 0; (unsigned)i < reatores.size(); i++){
 			distancia = Distancia(ponto,reatores[i].loc);
 			reatores[i].distance = distancia;
 			if(primeiro){
 				distanciamaisperto = distancia;
 				//strcpy(nomemaisperto,(*it).name);
 				primeiro = false;
 			}
 			else if(distancia < distanciamaisperto){
 				distanciamaisperto = distancia;
 				nomemaisperto = reatores[i].name;
 				//strcpy(nomemaisperto,(*it).name);
 			}
 		}
 		for(long double pot = PMIN; pot < PMAX && flux > FMIN; pot += PASSO){
 			x = pow(10,pot);
 			flux = 0;
 			for(int i = 0; (unsigned)i < reatores.size(); i++){
 				if(!reatores[i].country.compare("JP")){
 					fluxparticular = reatores[i].Fluxo(x);
 					if(oscila)
 						flux += survivalProb(0,0,reatores[i].distance * pow(10,-5)/*cm to km*/,x * pow(10,-3)/*MeV to GeV*/,true) * fluxparticular;
 					else
 						flux += fluxparticular;
 				}
 			}
 			saida << x << " " << flux << endl;
 		}
 		cout << "File " << nomesaida << " with the total flux on the coordinates " << ponto.lat << " " << ponto.lon << " successfully created.";
 		cout << " The nearest reactor to this location is " << nomemaisperto << " (" << distanciamaisperto/100 << "m)." << endl; 
 		cout << "Do you want to colect more data? (1=y, 0=n): ";
 		cin >> mais;
 		saida.close();
 		saida.clear();
 	}while(mais);
 	return(0);
}
