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
	double Pee = survivalProb_3genSun_Ad_freeparam(Enu, Rfrac, osc_p);

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
		for(int j = 0; j < detec.numberOfBins; j++)
			N[i*detec.numberOfBins + j] = op.exposure[i]*N_1d[j];
	}
}

int main(){
	/***********************/	
	/*OUTPUT FILE DIRECTORY*/
	/***********************/
	string outdir = "results/";

	ofstream saida;
	string file_name_base = "junoevents_regions_oscillation_daynight_";

	int n_E_Bins = 140;
	int n_cosz_Bins = 1;
	double binMin = 2;
	double binMax = 16;
	double Eres = 3;

	double *T = new double [n_E_Bins*n_cosz_Bins];
	double *N_obs = new double [n_E_Bins*n_cosz_Bins];
	double *N_BG = new double [n_E_Bins*n_cosz_Bins];
	double *N_BG_small = new double [n_E_Bins];
	double *exposure = new double [n_cosz_Bins];

	zprime_param zp;
	zp.antinu = false;
	zp.gX = zp.a = zp.e = zp.mZprime = zp.QlL = zp.QlR = zp.QnuL = 0;

	string exposure_file = "JUNO_exposure_100bins.txt";
	seta_exposure(exposure, n_cosz_Bins, exposure_file);
	seta_dados_Ne();

	osc_param osc_p;
	osc_p.seta_osc_param(false);
	osc_p.n_exp_bins = n_cosz_Bins;
	osc_p.exposure = exposure;

	osc_p.theta12 = asin(sqrt(0.307));
	osc_p.Dmq21 = 4.8e-5;
	osc_p.constroi_PMNS();
	
	signal_param sig;
	sig.flux = &Flux;
	sig.max_flux_E = 18.79; /*maximum solar flux energy*/
	sig.xsection = &Xsection;
	sig.xsection_param = (void *)(&zp);
	sig.oscillation_prob = &survivalProb_param;
	sig.oscillation_param = (void *)(&osc_p);

	detector_param detec;
	detec.numberOfBins = n_E_Bins;
	detec.Ebinmin = binMin;
	detec.Ebinmax = binMax;
	detec.binsize = (binMax - binMin)/n_E_Bins;
	detec.Eres = Eres;
	detec.Ntau = &Ntau_func;

	calcula_eventos_Bins_duplos(T, sig, detec);
	string BG_file = "JUNO_BG_total.txt";
	seta_BG(N_BG_small, n_E_Bins, BG_file);
	for(int i = 0; i < n_E_Bins; i++){
		for(int j = 0; j < n_cosz_Bins; j++){
			N_BG[i + j*n_E_Bins] = (1.0/n_cosz_Bins)*N_BG_small[i];
		}
	}

	stat_params p;
	p.nBins = n_E_Bins*n_cosz_Bins;
	p.T = T;
	p.BG = N_BG;

	p.sd = 0.05;
	p.sb = 0.15;

	/////////////////

	double siq_ini = 0.2, siq_fim = 0.45;
	double Dmq_ini = 0, Dmq_fim = 1.2e-4;
	
	int n_siq = 50, n_Dmq = 60;

	int n_of_arqs, n_arq_atual;
	//cout << "Entre com o número de arquivos: ";
	cin >> n_of_arqs;
	//cout << "Entre com o número do arquivo atual: ";
	cin >> n_arq_atual;

	stringstream ss;
	string numero_arquivo;
	string txt = ".txt";
	ss << n_arq_atual;
	ss >> numero_arquivo;

	saida.open(outdir + file_name_base + numero_arquivo + txt);
	
	double siq_step = (siq_fim - siq_ini)/n_siq;
	double Dmq_step = (Dmq_fim - Dmq_ini)/n_Dmq;

	int n_siq_piece = n_siq/n_of_arqs;
	double siq_ini_desse = siq_ini + (n_arq_atual-1)*n_siq_piece*siq_step;
	double siq_fim_desse = siq_ini + n_arq_atual*n_siq_piece*siq_step;

	for(double siq = siq_ini_desse; siq < siq_fim_desse; siq += siq_step){
		for(double Dmq = Dmq_ini; Dmq < Dmq_fim; Dmq += Dmq_step){
			osc_p.theta12 = asin(sqrt(siq));
			osc_p.Dmq21 = Dmq;
			osc_p.constroi_PMNS();
			sig.oscillation_param = (void *)(&osc_p);
			calcula_eventos_Bins_duplos(N_obs, sig, detec);
			for(int i = 0; i < n_E_Bins*n_cosz_Bins; i++)
				N_obs[i] += N_BG[i];
			p.N_obs = N_obs;
			saida << siq << " " << Dmq << " " << chisquare_twopulls_JUNO(p, false) << endl;
		}
		saida << endl;
	}	

	delete [] T;
	delete [] N_obs;
	delete [] N_BG;
	delete [] N_BG_small;
	delete [] exposure;

	saida.close();
	return(0);
}