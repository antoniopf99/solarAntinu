#include<iostream>
#include<cmath>
#include<fstream>
#include"solarflux.h"
#include"crosssections.h"
#include"detectors.h"
#include"general_adiabatic.h"
#include"constant_matter_evolution.h"
#include"statistics.h"

using namespace std;

double Xsection(double Enu, double T, void *p){
	zprime_param zp = *(zprime_param *)p;
	zp.Enu = Enu;
	zp.T = T;
	return(nueScatteringBSM(zp));
}

double survivalProb(double Enu, void *p){
	osc_param op = *(osc_param *)(p);
	if(op.prob_ee)
		return(sun_earth_prob_2layers(0.05, Enu, op));
	return(1 - sun_earth_prob_2layers(0.05, Enu, op));
}

double Flux(double Enu){
	return(solarFlux(3, Enu));
}

double countRate_SM(double T, double cosz){
	signal_param sig;
	
	sig.flux = &Flux;
	sig.max_flux_E = 18.79;
	
	sig.xsection = &Xsection;
	sig.oscillation_prob = &survivalProb;

	zprime_param xp;
	xp.antinu = false;
	xp.gX = xp.a = xp.e = xp.mZprime = xp.QlL = xp.QlR = xp.QnuL = 0;

	sig.xsection_param = (void *)(&xp);
	
	osc_param osc_p;
	osc_p.cosz = cosz;
	osc_p.seta_osc_param(false);
	osc_p.constroi_PMNS();

	sig.oscillation_param = (void *)(&osc_p);

	return(countRate_nu_e_Scattering(T, sig));
}

double Ntau_func(double E){
	//double Nelectrons = 6*(1e23);
	double Nelectrons = 1*3.38*pow(10,32);
	double year = 3.154*pow(10,7);
	year *= 1;
	double Ntau = Nelectrons*year;
	return(Ntau);
}

int main(){
	/***********************/	
	/*OUTPUT FILE DIRECTORY*/
	/***********************/
	string outdir = "results/";
	ofstream saida;
	string file = "compare_sginal.txt";
	//saida.open(outdir + file);

	seta_dados_Ne();

	for(double E = 0.01; E < 20; E += 0.1)
		for(double cosz = -1; cosz < 1; cosz += 0.01)
			cout << E << " " << cosz << " " << Ntau_func(0)*countRate_SM(E, cosz) << endl;

	saida.close();
	return(0);
}