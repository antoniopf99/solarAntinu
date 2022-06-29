#include"detectors.h"

double energyResolution(double Eobs, double Ereal, double Res){
	if(Ereal <= 0){
		cout << "ERROR COMPUTING EVENT AT ZERO ENERGY! (At the energyResolution function)." << endl;
		return(0);
	}
	double d = ( (Res/100) * sqrt(Ereal) );// / (2 * sqrt(2*log(2)));
	double N = 1 / (sqrt(2*M_PI) * d);
	return(N * exp(-0.5 * pow((Eobs - Ereal)/d,2)));
}

double resolutionFactor(double Eini, double Efinal, double Ereal, double Res){
	double d = ( (Res/100) * sqrt(Ereal) );
	double final = gsl_sf_erf((Efinal - Ereal)/(sqrt(2)*d));
	double ini = gsl_sf_erf((Eini - Ereal)/(sqrt(2)*d));
	return((final - ini)/(2*sqrt(2*M_PI)));
}

double countRate_nu_e_Scattering_integrand(double E, void *sig_param){
	signal_param p = *(signal_param *)(sig_param);
	double T = ((zprime_param *)p.xsection_param)->T;
	double flux = p.flux(E);
	double xsection = p.xsection(E, T, p.xsection_param);
	double survival_prob = p.oscillation_prob(E, p.oscillation_param);
	osc_param op = *(osc_param*)(p.oscillation_param);
	if(isnan(survival_prob)){
		cout << "SURVIVAL PROB TÁ NAN! EM:" << endl;
		cout << "E = " << E << endl;
		cout << "cosz = " << op.cosz << endl;
		cout << endl;
	}
	return(flux*xsection*survival_prob);
}

double countRate_nu_e_Scattering(double T, signal_param sig){
	double Eanalmin, Eanalmax;
	Eanalmin = 0.5*(T + sqrt(T*(T + 2*0.511)));
	Eanalmax = sig.max_flux_E;

	zprime_param zp = *(zprime_param *)(sig.xsection_param);
	osc_param osc_p = *(osc_param *)(sig.oscillation_param);

	zp.T = T;
	sig.xsection_param = (void *)(&zp);

	integration_parameters quad_param;
	quad_param.sig = sig;

	gsl_function F;
	F.function = &countRate_nu_e_Scattering_integrand;

	double integral, resultado = 0, abserr;
	size_t neval;

	//gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

	/*contribution from the electron neutrino flux*/
	zp.NC_and_CC = true;
	osc_p.prob_ee = true;
	quad_param.sig.xsection_param = (void *)(&zp);
	quad_param.sig.oscillation_param = (void *)(&osc_p);
	F.params = (void *)(&quad_param);
	//gsl_integration_qags(&F, Eanalmin, Eanalmax, ERRO_ABS, ERRO_REL, MAX_SUB_INT, w, &integral, &abserr);
	gsl_integration_qng(&F, Eanalmin, Eanalmax, ERRO_ABS, ERRO_REL, &integral, &abserr, &neval);
	resultado += integral;

	/*contribution from the heavy neutrino flavors flux*/
	zp.NC_and_CC = false;
	osc_p.prob_ee = false;
	quad_param.sig.xsection_param = (void *)(&zp);
	quad_param.sig.oscillation_param = (void *)(&osc_p);
	F.params = (void *)(&quad_param);
	//gsl_integration_qags(&F, Eanalmin, Eanalmax, ERRO_ABS, ERRO_REL, MAX_SUB_INT, w, &integral, &abserr);
	gsl_integration_qng(&F, Eanalmin, Eanalmax, ERRO_ABS, ERRO_REL, &integral, &abserr, &neval);
	resultado += integral;

	//gsl_integration_workspace_free(w);
	return(resultado);
}

double calcula_eventos_Bins_integrand(double T, void *quad_param){
	integration_parameters qp = *(integration_parameters *)(quad_param);

	double h = ((qp.Tanalmax - qp.Tanalmin) / (qp.N-1));

	int i = (int)((T - qp.Tanalmin)/h);
	double x1 = qp.T[i];
	double x2 = qp.T[i+1];
	double y1 = qp.countRate[i];
	double y2 = qp.countRate[i+1];

	double countRate = ((y2-y1)/(x2-x1))*(T-x1) + y1;

	return(countRate*resolutionFactor(qp.E_current_bin, qp.E_current_bin + qp.detec.binsize, T, qp.detec.Eres));
}

void calcula_eventos_Bins(double *Bins, signal_param sig, detector_param detec){

	/*NUMBER OF POINTS IN COUNT RATE INTERPOLATION*/
	int N = 1000;

	/*Analysis energy range of THE TRUE electron recoil*/
	double Tanalmin = 0.001;
	double Tanalmax = 20;
	double h = (Tanalmax - Tanalmin) / N;

	double *countRate = new double [N+1];

	for(int i = 0; i <= N; i++){
		//printf("Computando, %.2f %\n", 100*(float)(i+1)/(N+1));
		countRate[i] = countRate_nu_e_Scattering(Tanalmin + i*h, sig);
	}

	double E = detec.Ebinmin;
	for(int i = 0; i < detec.numberOfBins; i++){
		Bins[i] = 0;
		/*Simpsons integration:*/
		/////
		Bins[i] += countRate[0]*resolutionFactor(E, E + detec.binsize, Tanalmin, detec.Eres);
		Bins[i] += countRate[N]*resolutionFactor(E, E + detec.binsize, Tanalmin + N*h, detec.Eres);
		for(int j = 0; j < N; j+=2){
			Bins[i] += 4*countRate[j]*resolutionFactor(E, E + detec.binsize, Tanalmin + j*h, detec.Eres);
			Bins[i] += 2*countRate[j+1]*resolutionFactor(E, E + detec.binsize, Tanalmin + (j+1)*h, detec.Eres);
		}
		Bins[i] *= h/3;
		/////
		E += detec.binsize;
		Bins[i] *= detec.Ntau(E);
	}
	delete [] countRate;
}

//void calcula_eventos_Bins(double *Bins, signal_param sig, detector_param detec){
//	/*NUMBER OF POINTS IN COUNT RATE INTERPOLATION*/
//	int N = 1000;
//	/*Analysis energy range of THE TRUE electron recoil*/
//	double Tanalmin = 0.001;
//	double Tanalmax = 20;
//	double h = (Tanalmax - Tanalmin) / N;
//	double *T = new double [N+1];
//	double *countRate = new double [N+1];
//	for(int i = 0; i <= N; i++){
//		//printf("Computando, %.2f %\n", 100*(float)(i+1)/(N+1));
//		T[i] = Tanalmin +i*h;
//		countRate[i] = countRate_nu_e_Scattering(T[i], sig);
//	}
//	integration_parameters quad_param;
//	quad_param.sig = sig;
//	quad_param.detec = detec;
//	quad_param.Tanalmin = Tanalmin;
//	quad_param.Tanalmax = Tanalmax;
//	quad_param.N = N+1;
//	quad_param.T = T;
//	quad_param.countRate = countRate;
//	gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
//	gsl_function F;
//	F.function = &calcula_eventos_Bins_integrand;
//	double integral, resultado = 0, abserr;
//	size_t neval;
//	double E = detec.Ebinmin;
//	//quad_param.E_current_bin = detec.Ebinmin;
//	for(int i = 0; i < detec.numberOfBins; i++){
//		Bins[i] = 0;
//		//F.params = &quad_param;
//		//gsl_integration_qng(&F, Tanalmin, Tanalmax, ERRO_ABS, ERRO_REL, &integral, &abserr, &neval);
//		//gsl_integration_qags(&F, Tanalmin, Tanalmax, ERRO_ABS, ERRO_REL, MAX_SUB_INT, w, &integral, &abserr);
//		//Bins[i] = integral;
//		/*Simpsons integration:*/
//		/////
//		Bins[i] += countRate[0]*resolutionFactor(E, E + detec.binsize, Tanalmin, detec.Eres);
//		Bins[i] += countRate[N]*resolutionFactor(E, E + detec.binsize, Tanalmin + N*h, detec.Eres);
//		for(int j = 0; j < N; j+=2){
//			Bins[i] += 4*countRate[j]*resolutionFactor(E, E + detec.binsize, Tanalmin + j*h, detec.Eres);
//			Bins[i] += 2*countRate[j+1]*resolutionFactor(E, E + detec.binsize, Tanalmin + (j+1)*h, detec.Eres);
//		}
//		Bins[i] *= h/3;
//		/////
//		E += detec.binsize;
//		quad_param.E_current_bin += detec.binsize;
//		Bins[i] *= detec.Ntau(quad_param.E_current_bin);
//	}
//	gsl_integration_workspace_free(w);
//	delete [] T;
//	delete [] countRate;
//}

void seta_BG(double *BG, int numberOfBins, string bg_file_name){
	/*background file directory*/
	string data_dir = "data/";
	ifstream bg_file;
	bg_file.open(data_dir + bg_file_name);
	for(int i = 0; i < numberOfBins; i++)
		bg_file >> BG[i];
	bg_file.close();
}

void seta_exposure(double *exposure, int numberOfExpBins, string exp_file_name){
	/*exposure file directory*/
	string data_dir = "data/";
	ifstream exp_file;
	exp_file.open(data_dir + exp_file_name);
	for(int i = 0; i < numberOfExpBins; i++)
		exp_file >> exposure[i];
	exp_file.close();
}