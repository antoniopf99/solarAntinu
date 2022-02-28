#include<iostream>
#include"solarmodel.h"

vector<data_point> Ne_data;

long double interpola_electron_density(long double Rfrac){
	long double x1, x2, y1, y2, resultado;
	bool foi = false;

	/****************/
	/*DATA FILE NAME*/
	/****************/
	string datafile = "data/SSMelectrondensity.txt";

	ifstream arquivo;
	arquivo.open(datafile);
	arquivo >> x1;
	arquivo >> y1;
	if(x1 > Rfrac){
 		resultado = y1;
 		foi = true;
 	}
	while(!arquivo.eof() && !foi){
 		arquivo >> x2;
 		arquivo >> y2;
 		if(x1 <= Rfrac && Rfrac <= x2){
 			resultado = ((y2-y1)/(x2-x1))*(Rfrac-x1) + y1;
 		 	foi = true;
 		}
 		x1 = x2;
 		y1 = y2;
 		arquivo.peek();
 	}
 	arquivo.close();
 	if(!foi)
 		resultado = y1;

 	return(resultado);
}

void seta_dados_Ne(){
	data_point data;

	for(double Rfrac = 0; Rfrac < 1; Rfrac += 0.001){
		data.Rfrac = Rfrac;
		data.dat = interpola_electron_density(Rfrac);
		Ne_data.push_back(data);
	}
}

long double electron_density_sun_NA(long double Rfrac){
	int size = Ne_data.size();
	
	if(Rfrac == 1)
		return((Ne_data[size - 1]).dat);
	
	int i = (int)(size * Rfrac);
	double x1 = ((double)i)/((double)size);
	double x2 = ((double)(i+1))/((double)size);
	double y1 = Ne_data[i].dat;
	double y2 = Ne_data[i+1].dat;
	//cout << size << endl;
	return(((y2-y1)/(x2-x1))*(Rfrac-x1) + y1);
}

long double interpola_dados_solares(long double Rfrac, int dado){
	long double x1, x2, y1[12], y2[12], resultado;
	bool foi = false;

	/****************/
	/*DATA FILE NAME*/
	/****************/
	string datafile = "data/SSM.txt";

	ifstream arquivo;
	arquivo.open(datafile);
	arquivo >> x1;
	for(int i = 0; i < 12; i++)
		arquivo >> y1[i];
	if(x1 > Rfrac){
 		resultado = y1[dado];
 		foi = true;
 	}
	while(!arquivo.eof() && !foi){
 		arquivo >> x2;
 		for(int i = 0; i < 12; i++)
			arquivo >> y2[i];
 		if(x1 <= Rfrac && Rfrac <= x2){
 			resultado = ((y2[dado]-y1[dado])/(x2-x1))*(Rfrac-x1) + y1[dado];
 		 	foi = true;
 		}
 		x1 = x2;
 		for(int i = 0; i < 12; i++)
 			y1[i] = y2[i];
 		arquivo.peek();
 	}
 	arquivo.close();
 	if(!foi)
 		resultado = y1[dado];

 	return(resultado);
}