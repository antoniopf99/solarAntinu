#include<iostream>
#include<cmath>
#include<fstream>
#include"solarflux.h"
#include"crosssections.h"
#include"detectors.h"
#include"general_adiabatic.h"
#include"statistics.h"

using namespace std;

double Xsection(double Enu, double T, void *p){
	zprime_param zp = *(zprime_param *)p;
	zp.Enu = Enu;
	zp.T = T;
	return(nueScatteringBSM(zp));
}

double survivalProb(double Enu, void *p){
	//return(survival_prob_sun_general_adiabatic(0.05, Enu, 1, 0, 0));
	if(*(bool *)p)
		return(survivalProb_3genSun_Ad(Enu, 0.05));
	return(1 - survivalProb_3genSun_Ad(Enu, 0.05));
}

double Flux(double Enu){
	return(solarFlux(3, Enu));
}

double Ntau_func(double E){
	double Nelectrons1 = 7.9*3.38*pow(10,32);
	double Nelectrons2 = 12.2*3.38*pow(10,32);
	double Nelectrons3 = 16.5*3.38*pow(10,32);
	double year = 3.154*pow(10,7);
	year *= 10;
	double Ntau1 = Nelectrons1*year;
	double Ntau2 = Nelectrons2*year;
	double Ntau3 = Nelectrons3*year;
	double Ntau;
	if(E <= 3)
		Ntau = Ntau1;
	else if(E <= 5)
		Ntau = Ntau2;
	else
		Ntau = Ntau3;
	return(Ntau);
}

int main(){
	/***********************/	
	/*OUTPUT FILE DIRECTORY*/
	/***********************/
	string outdir = "results/";

	ofstream saida;
	string file = "junoevents_regions.txt";
	saida.open(outdir + file);

	int numberOfBins = 140;
	double binMin = 2;
	double binMax = 16;
	double N_SM[210];
	double N_obs[210];
	double N_BG[210];
	double Eres = 3;

	//double binsize = (binMax-binMin)/numberOfBins;

	zprime_param zp;
	zp.antinu = false;
	zp.gX = zp.a = zp.e = zp.mZprime = zp.QlL = zp.QlR = zp.QnuL = 0;

	bool prob_ee = true;

	seta_dados_Ne();

	signal_param sig;
	sig.flux = &Flux;
	sig.max_flux_E = 18.79; /*maximum solar flux energy*/
	sig.xsection = &Xsection;
	sig.xsection_param = (void *)(&zp);
	sig.oscillation_prob = &survivalProb;
	sig.oscillation_param = (void *)(prob_ee);

	detector_param detec;
	detec.numberOfBins = 140;
	detec.Ebinmin = binMin;
	detec.Ebinmax = binMax;
	detec.binsize = (binMax - binMin)/numberOfBins;
	detec.Eres = Eres;
	detec.Ntau = &Ntau_func;

	calcula_eventos_Bins(N_SM, sig, detec);
	string BG_file = "JUNO_BG_total.txt";
	seta_BG(N_BG, numberOfBins, BG_file);

	stat_params p;
	p.nBins = 140;
	p.T = N_SM;
	p.BG = N_BG;

	p.sd = 0.05;
	p.sb = 0.15;

	double Mz, e;
	for(double logMz = -2; logMz < 4; logMz += 0.1){
		for(double loge = -6; loge < -1; loge += 0.1){
			Mz = pow(10, logMz);
			e = pow(10, loge);
			zp.mZprime = Mz;
			zp.e = e;
			sig.xsection_param = (void *)(&zp);
			calcula_eventos_Bins(N_obs, sig, detec);
			for(int i = 0; i < numberOfBins; i++)
				N_obs[i] += N_BG[i];
			p.N_obs = N_obs;
			saida << Mz << " " << e << " " << chisquare_twopulls_JUNO(p, false) << endl;
		}
		saida << endl;
	}	

	saida.close();
	return(0);
}