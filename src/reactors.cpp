#include<iostream>
#include<cmath>
#include<fstream>
#include"reactors.h"

using namespace std;

long double Lambda(int i, long double E){

	/**********************************************/
	/*FLUX CALCULATION PARAMETERS: LAMBDA FUNCTION*/
	/**********************************************/
	long double a[4][6];
	/*235 U*/
	a[0][0] = 3.217;
	a[0][1] = -3.111;
	a[0][2] = 1.395;
	a[0][3] = -3.690 * pow(10,-1);
	a[0][4] = 4.445 * pow(10,-2);
	a[0][5] = -2.053 * pow(10,-3);

	/*238 U*/
	a[1][0] = 4.833 * pow(10,-1);
	a[1][1] = 1.927 * pow(10,-1);
	a[1][2] = -1.283 * pow(10,-1);
	a[1][3] = -6.762 * pow(10,-3);
	a[1][4] = 2.233 * pow(10,-3);
	a[1][5] = -1.536 * pow(10,-4);
	
	/*239 Pu*/
	a[2][0] = 6.413;
	a[2][1] = -7.432;
	a[2][2] = 3.535;
	a[2][3] = -8.820 * pow(10,-1);
	a[2][4] = 1.025 * pow(10,-1);
	a[2][5] = -4.550 * pow(10,-3);
	
	/*241 Pu*/
	a[3][0] = 3.251;
	a[3][1] = -3.204;
	a[3][2] = 1.428;
	a[3][3] = -3.675 * pow(10,-1);
	a[3][4] = 4.254 * pow(10,-2);
	a[3][5] = -1.896 * pow(10,-3);

	long double soma = 0;
	for(int j = 0; j < 6; j++)
		soma += (a[i][j] * pow(E,j));
	return(exp(soma));
}

long double FatorA(int i, long double E){

	/**************************************************************/
	/*FLUX CALCULATION PARAMETERS: ISOTOPE ENERGY, POWER FRACTIONS*/
	/**************************************************************/
	/*ISOTOPE ENERGY IN MeV*/
	long double Q[4];
	/*235 U*/
	Q[0] = 202.36;
	/*238 U*/
	Q[1] = 205.99;
	/*239 Pu*/
	Q[2] = 211.12;
	/*235 Pu*/
	Q[3] = 214.26;
	/*POWER FRACTIONS*/
	long double p[4][3];
	/*235 U*/
	p[0][0] = 0.561; /*PWR*/
	p[0][1] = 0; /*MOX*/
	p[0][2] = 0.543; /*PHWR*/
	/*238 U*/
	p[1][0] = 0.076; /*PWR*/
	p[1][1] = 0.081; /*MOX*/
	p[1][2] = 0.411; /*PHWR*/
	/*239 Pu*/
	p[2][0] = 0.307; /*PWR*/
	p[2][1] = 0.708; /*MOX*/
	p[2][2] = 0.022; /*PHWR*/
	/*241 Pu*/
	p[3][0] = 0.056; /*PWR*/
	p[3][1] = 0.212; /*MOX*/
	p[3][2] = 0.024; /*PHWR*/

	long double soma = 0;
	for(int j = 0; j < 4; j++)
		soma += (p[j][i] / /*emKeV*/Q[j]) * ((Lambda(j,/*emMeV*/E))/* / 1000*/) ; // converter 1/MeV para 1/KeV!!
	return(soma);
}

void leOsDados(vector<Reactor> &reatores, int quantosAnos){

	/*********************/
	/*DATA DIRECTORY NAME*/
	/*********************/
	string datadir = "data/";

	ifstream arquivo;
	string nomearq;
	Reactor reat;
	double lf, fraction;
	bool tem;
	do{
		cout << "Enter with the name of the file with reactor data: ";
		cin >> nomearq;
		arquivo.open(datadir + nomearq);
		if(!arquivo.is_open())
			cout << "Could not open that file, try again." << endl;
 	}while(!arquivo.is_open());
 	cout << "Enter with the fraction of the year that the detector was operational: ";
 	cin >> fraction;
	while(!arquivo.eof()){
 		arquivo >> reat.country >> reat.name;
 		arquivo >> reat.loc.lat >> reat.loc.lon;
 		arquivo >> reat.tipo;
 		if(!reat.tipo.compare("PHWR"))
 			reat.type = 2;
 		else
 			reat.type = 0;
 		arquivo >> reat.MOX;
 		arquivo >> reat.thermalPower;
 		reat.loadFactor = 0;
 		for(int i = 0; i < 12; i++){
 			arquivo >> lf;
 			reat.loadFactor += lf / 12;
 		}
 		reat.loadFactor *= fraction;
 		tem = false;
 		for(int i = 0; i < (int)reatores.size() && !tem; i++){
 			if(!reat.country.compare(reatores[i].country) && !reat.name.compare(reatores[i].name)){
 				reatores[i].loadFactor += reat.loadFactor;
 				tem = true;
 			}
 		}
 		if(!tem)
 			reatores.push_back(reat);
 		arquivo.peek();
 	}
 	arquivo.close();
}

long double customFlux(int type, int MOX, long double thermalPower, long double loadFactor, long double E, long double distance){
	long double base = ((/*emMeVs*/(6.2415093433 * pow(10,18) * thermalPower) * (loadFactor/100) * FatorA(type,E)) / (4 * M_PI * pow(distance,2)));
	if(!MOX)
		return(base);
	long double mox = ((/*emMeVs*/(6.2415093433 * pow(10,18) * thermalPower) * (loadFactor/100) * FatorA(1,E)) / (4 * M_PI * pow(distance,2)));
	return(0.7 * base + 0.3 * mox);
}

long double Reactor::Fluxo(long double E){
	return(customFlux(type, MOX, thermalPower, loadFactor, E, distance));
}