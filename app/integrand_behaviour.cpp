#include<iostream>
#include<cmath>
#include<fstream>
#include<sstream>
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

double survivalProb_param(double Enu, void *p){
	osc_param osc_p = *(osc_param *)p;
	double Rfrac = 0.05;

	//double Pee = osc_p.exposure[osc_p.i]*survivalProb_3genSun_Ad_freeparam(Enu, Rfrac, osc_p);//*sun_earth_prob_2layers(Rfrac, Enu, osc_p);
	//double Pee = survivalProb_3genSun_Ad_freeparam(Enu, Rfrac, osc_p);
	double Pee = sun_earth_prob_5layers(Rfrac, Enu, osc_p);

	if(osc_p.prob_ee)
		return(Pee);
	return(1 - Pee);
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

void calcula_eventos_Bins_duplos(double *N, signal_param sig, detector_param detec){
	osc_param op = *(osc_param *)sig.oscillation_param;
	double binsize = 2.0/op.n_exp_bins;
	//op.cosz = binsize/2 - 1;
	op.cosz = -1;
	double N_1d[210];
	for(int i = 0; i < op.n_exp_bins; i++){
		op.i = i;
		sig.oscillation_param = (void *)(&op);
		calcula_eventos_Bins(N_1d, sig, detec);
		for(int j = 0; j < detec.numberOfBins; j++){
			N[i*detec.numberOfBins + j] = op.exposure[i]*N_1d[j];
			//N[i*detec.numberOfBins + j] = (1.0/op.n_exp_bins)*N_1d[j];
		}
		op.cosz += binsize;
	}
}

int main(){
	seta_dados_Ne();
	/***********************/	
	/*OUTPUT FILE DIRECTORY*/
	/***********************/
	string outdir = "results/";

	ofstream saida;
	string file_name = "junoevents_daynight.txt";
	saida.open(outdir + file_name);

	int n_E_Bins = 140;
	int n_cosz_Bins = 100;
	double binMin = 2;
	double binMax = 16;
	double binsize = (binMax - binMin)/n_E_Bins;
	double Eres = 3;

	double Tanalmin = 0.001;
	double Tanalmax = 20;

	zprime_param zp;
	zp.antinu = false;
	zp.gX = zp.a = zp.e = zp.mZprime = zp.QlL = zp.QlR = zp.QnuL = 0;

	osc_param osc_p;
	osc_p.seta_osc_param(false);

	osc_p.theta12 = 0.587252;//asin(sqrt(0.307));
	osc_p.Dmq21 = 4.8e-5;
	osc_p.cosz = 0.36;
	osc_p.constroi_PMNS();

	detector_param detec;
	detec.numberOfBins = n_E_Bins;
	detec.Ebinmin = binMin;
	detec.Ebinmax = binMax;
	detec.binsize = binsize;
	detec.Eres = Eres;
	detec.Ntau = &Ntau_func;
	
	signal_param sig;
	sig.flux = &Flux;
	sig.max_flux_E = 18.79; /*maximum solar flux energy*/
	sig.xsection = &Xsection;
	sig.xsection_param = (void *)(&zp);
	sig.oscillation_prob = &survivalProb_param;
	sig.oscillation_param = (void *)(&osc_p);

	int N = 1000;
	double h = (Tanalmax - Tanalmin) / N;

	double *T = new double [N+1];	
	double *countRate = new double [N+1];
	for(int i = 0; i <= N; i++){
		T[i] = Tanalmin + i*h;
		countRate[i] = countRate_nu_e_Scattering(T[i], sig);
	}

	integration_parameters quad_param;
	quad_param.sig = sig;
	quad_param.detec = detec;
	quad_param.Tanalmin = Tanalmin;
	quad_param.Tanalmax = Tanalmax;
	quad_param.N = N;
	quad_param.T = T;
	quad_param.countRate = countRate;

	for(double Trec = Tanalmin; Trec < Tanalmax; Trec += 0.01)
		cout << Trec << " " << countRate_nu_e_Scattering(Trec, sig) << " " << calcula_eventos_Bins_integrand(Trec, (void*)(&quad_param)) << endl;

	//for(double T = Tanalmin; T < Tanalmax; T += 0.1){
	//	for(double E = 0.5*(T + sqrt(T*(T + 2*0.511))); E < sig.max_flux_E; E += 0.1){
	//		zp.T = T;
	//		sig.xsection_param = (void *)(&zp);
	//		cout << T << " " << E << " " << countRate_nu_e_Scattering_integrand(E, (void*)(&sig)) << endl;
	//	}
	//}

	delete [] T;
	delete [] countRate;
	saida.close();
	return(0);
}