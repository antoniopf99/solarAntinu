#include<iostream>
#include<cmath>
#include<fstream>
#include"magnus_evolution.h"
#include"constant_matter_evolution.h"
#include"detectors.h"

using namespace std;

double survivalProb(double Enu, void *p){
	osc_param op = *(osc_param *)(p);
	double Pee = 0;
	//double binsize = 2.0/op.n_exp_bins;
	//double cosz = binsize/2;
	for(int i = 0; i < op.n_exp_bins; i++){
		//Pee += op.exposure[i]*survivalProb_3genSun_Ad_freeparam(Enu, Rfrac, op);//*sun_earth_prob_2layers(0.05, cosz, Enu, op);
		//cosz += binsize;
	}
	return(Pee);
}



int main(){
	seta_dados_Ne();
	
	string exposure_file = "JUNO_exposure_100bins.txt";
	double exposure[100];

	osc_param op;
	seta_exposure(exposure, 100, exposure_file);
	op.n_exp_bins = 100;
	op.exposure = exposure;
	op.seta_osc_param(false);
	op.delta_CP = -M_PI/2;
	op.constroi_PMNS();

	op.h = 1;

	double Rfrac = 0.05;

	//double E;
	//for(double logE = -2; logE < 8; logE += 0.05){
	//	E = pow(10, logE);
	//	cout << E << " " << survival_prob_sun_magnus(Rfrac, E, op) << endl;
	//}

	//for(double E = 200; E < 1800; E += 0.5){
	//	cout << E << " ";
	//	op.cosz = cos(M_PI/10.0);
	//	cout << atmospheric_neutrino_5layers(1, 0, E, op) << " ";
	//	op.cosz = cos(M_PI/4.0);
	//	cout << atmospheric_neutrino_5layers(1, 0, E, op) << " ";
	//	op.cosz = cos(M_PI/3.0);
	//	cout << atmospheric_neutrino_5layers(1, 0, E, op) << " ";
	//	op.cosz = cos(M_PI/2.5);
	//	cout << atmospheric_neutrino_5layers(1, 0, E, op) << endl;
	//}

	double probmagnus, probconst5;
	for(double E = 0.1; E < 20; E += 0.01){
			cout << E << " ";
			op.cosz = 0.1;
			probmagnus = survival_prob_sun_ad_earth_magnus(0.05, E, op);
			probconst5 = sun_earth_prob_5layers(Rfrac, E, op);
			cout << probmagnus << " ";
			cout << probconst5 << " ";
			cout << probconst5 - probmagnus << " ";
			op.cosz = 0.4;
			probmagnus = survival_prob_sun_ad_earth_magnus(0.05, E, op);
			probconst5 = sun_earth_prob_5layers(Rfrac, E, op);
			cout << probmagnus << " ";
			cout << probconst5 << " ";
			cout << probconst5 - probmagnus << " ";
			op.cosz = 0.6;
			probmagnus = survival_prob_sun_ad_earth_magnus(0.05, E, op);
			probconst5 = sun_earth_prob_5layers(Rfrac, E, op);
			cout << probmagnus << " ";
			cout << probconst5 << " ";
			cout << probconst5 - probmagnus << " ";
			op.cosz = 0.9;
			probmagnus = survival_prob_sun_ad_earth_magnus(0.05, E, op);
			probconst5 = sun_earth_prob_5layers(Rfrac, E, op);
			cout << probmagnus << " ";
			cout << probconst5 << " ";
			cout << probconst5 - probmagnus << " ";
			op.cosz = 1;
			probmagnus = survival_prob_sun_ad_earth_magnus(0.05, E, op);
			probconst5 = sun_earth_prob_5layers(Rfrac, E, op);
			cout << probmagnus << " ";
			cout << probconst5 << " ";
			cout << probconst5 - probmagnus << endl;
	}

	//double E, prob;
	//for(double cosz = 0.7; cosz < 1; cosz += 0.001){
	//	for(double logE = -1; logE < 1.3; logE += 0.001){
	//		op.cosz = cosz;
	//		E = 1000*pow(10, logE);
	//		prob = atmospheric_neutrino_5layers(1, 1, E, op);
	//		//if(prob > 1 || prob < 0){
	//		//	cout << "PROB " << prob << " EM cosz = " << cosz << " e E = " << E/1000.0 << " GeV." << endl;
	//		//}
	//		cout << (E/1000.0) << " " << cosz << " " << prob << endl;
	//	}
	//	//cout << endl;
	//}


	return(0);
}
