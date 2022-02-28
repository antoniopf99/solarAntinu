#include<iostream>
#include<cmath>
#include<fstream>
#include<gsl/gsl_integration.h>
#include"crosssections.h"
#include"reactors.h"
#include"oscilation.h"
#include"detectors.h"

using namespace std;

double Xsection(double Enu, void *parametros){

	double T = ((double *)parametros)[0];
	double gBL = ((double *)parametros)[1];
	double mZprime = ((double *)parametros)[2];
	
	return(antinueScatteringBSM(1,gBL,0,0,mZprime,-1,-1,-1,Enu,T));
}

double normalizedXsection(double gBL, double mZprime, double T){
	double parametros[3] ={T, gBL, mZprime};

	gsl_function F;
	F.function = &Xsection;
	F.params = parametros;

	gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

	double result, error;

	double me = 0.510998; /*electron mass in MeV*/
	double Enumin = 0.5*T +sqrt(T*(T+2*me));

	gsl_integration_qag(&F, Enumin, 100, 0, 0.000001, 1000, 6, w, &result, &error);
	gsl_integration_workspace_free(w);

	return(result);
}

int main(){

	/***********************/	
	/*OUTPUT FILE DIRECTORY*/
	/***********************/
	string outdir = "results/";
	ofstream saida;
	string file = "XsectionSM.txt";
	saida.open(outdir + file);

	for(double T = 0; T < 12; T+=0.001)
		saida << T << " " << antinueScatteringSM(10, T) << endl;

//	double T;
//	double gBL = 5*pow(10,-4);
//
//	for(double exp = expmin; exp <= expmax; exp += step){
//		T = pow(10,exp);
//		saida << T << " " << pow(10,-45)*normalizedXsection(0,0,T) << endl;
//	}
//
//	saida.close();
//
//	file = "XsectionB-L100.txt";
//	saida.open(outdir + file);
//
//	for(double exp = expmin; exp <= expmax; exp += step){
//		T = pow(10,exp);
//		saida << T << " " << pow(10,-45)*normalizedXsection(gBL,100,T) << endl;
//	}
//
//	saida.close();
//
//	file = "XsectionB-L10.txt";
//	saida.open(outdir + file);
//
//	for(double exp = expmin; exp <= expmax; exp += step){
//		T = pow(10,exp);
//		saida << T << " " << pow(10,-45)*normalizedXsection(gBL,10,T) << endl;
//	}
//
//	saida.close();
//
//	file = "XsectionB-Lp1.txt";
//	saida.open(outdir + file);
//
//	for(double exp = expmin; exp <= expmax; exp += step){
//		T = pow(10,exp);
//		saida << T << " " << pow(10,-45)*normalizedXsection(gBL,0.1,T) << endl;
//	}
//
	saida.close();

	return(0);
}