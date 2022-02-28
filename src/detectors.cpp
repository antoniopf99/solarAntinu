#include"detectors.h"

long double energyResolution(long double Eobs, long double Ereal, double Res){
	if(Ereal <= 0){
		cout << "ERROR COMPUTING EVENT AT ZERO ENERGY! (At the energyResolution function)." << endl;
		return(0);
	}
	long double d = ( (Res/100) * sqrt(Ereal) );// / (2 * sqrt(2*log(2)));
	long double N = 1 / (sqrt(2*M_PI) * d);
	return(N * exp(-0.5 * pow((Eobs - Ereal)/d,2)));
}

double resolutionFactor(double Eini, double Efinal, double Ereal, double Res){
	double d = ( (Res/100) * sqrt(Ereal) );
	double final = gsl_sf_erf((Efinal - Ereal)/(sqrt(2)*d));
	double ini = gsl_sf_erf((Eini - Ereal)/(sqrt(2)*d));
	return((final - ini)/(2*sqrt(2*M_PI)));
	//int N = 500;
	//long double h = (Efinal - Eini) / N;
	//long double soma = energyResolution(Eini,Ereal,Res) + energyResolution(Eini + N*h,Ereal,Res);
	//for(long i = 1; i < N; i+=2){
	//	soma += 4*energyResolution(Eini + i*h,Ereal,Res) + 2*energyResolution(Eini + (i+1)*h,Ereal,Res);
	//}
	//soma *= h/3;
	//return(soma);
}

long double countRate_nu_e_Scattering(
	long double (*Flux)(long double), 
	long double (*Xsection)(long double, long double, double *), 
	double *XsectionParameters, 
	long double (*survivalProb)(long double, double *), 
	double *survivalProbParameters,
	double T)
{
	long double soma, Eanalmin, Eanalmax, h;
	int N = 1000;
	Eanalmin = 0.5*(T + sqrt(T*(T + 2*0.511)));
	Eanalmax = 18.79;
	h = (Eanalmax - Eanalmin) / N;
	if(h<=0)
		return(0);
	soma = 0;
	soma += Flux(Eanalmin)*Xsection(Eanalmin, T, XsectionParameters)*survivalProb(Eanalmin, survivalProbParameters);
	soma += Flux(Eanalmin + N*h)*Xsection(Eanalmin + N*h, T, XsectionParameters)*survivalProb(Eanalmin + N*h, survivalProbParameters);
	for(int i = 1; i < N; i+=2){
		soma += 4*Flux(Eanalmin + i*h)*Xsection(Eanalmin + i*h, T, XsectionParameters)*survivalProb(Eanalmin + i*h, survivalProbParameters);
		soma += 2*Flux(Eanalmin + (i+1)*h)*Xsection(Eanalmin + (i+1)*h, T, XsectionParameters)*survivalProb(Eanalmin + (i+1)*h, survivalProbParameters);
	}
	soma *= h/3;
	return(soma);
}

double events_in_bin(
	long double (*Flux)(long double), 
	long double (*Xsection)(long double, long double, double *), 
	double *XsectionParameters, 
	long double (*survivalProb)(long double, double *), 
	double *survivalProbParameters,
	double Ebinmin,
	double binsize,
	long double Ntau,
	double Eres)
{
	long double soma, Eanalmin, Eanalmax, h;
	int N = 500;
	Eanalmin = 0.001;
	Eanalmax = 20;
	h = (Eanalmax - Eanalmin) / N;
	soma = 0;
	soma += countRate_nu_e_Scattering(Flux, Xsection, XsectionParameters, survivalProb, survivalProbParameters, Eanalmin)*resolutionFactor(Ebinmin, Ebinmin + binsize, Eanalmin, Eres);
	soma += countRate_nu_e_Scattering(Flux, Xsection, XsectionParameters, survivalProb, survivalProbParameters, Eanalmin + N*h)*resolutionFactor(Ebinmin, Ebinmin + binsize, Eanalmin + N*h, Eres);
	for(int i = 0; i < N; i+=2){
		soma += 4*countRate_nu_e_Scattering(Flux, Xsection, XsectionParameters, survivalProb, survivalProbParameters, Eanalmin + i*h)*resolutionFactor(Ebinmin, Ebinmin + binsize, Eanalmin + i*h, Eres);
		soma += 2*countRate_nu_e_Scattering(Flux, Xsection, XsectionParameters, survivalProb, survivalProbParameters, Eanalmin + (i+1)*h)*resolutionFactor(Ebinmin, Ebinmin + binsize, Eanalmin + (i+1)*h, Eres);
	}
	soma *= h/3;
	soma *= Ntau;
	return(soma);
}

void calcula_eventos_Bins(
	double *Bins,
	long double Ebinmin,
	long double Ebinmax,
	int numberOfBins,
	long double (*Ntau)(long double),
	long double (*Flux)(long double), 
	long double (*Xsection)(long double, long double, double *), 
	double *XsectionParameters, 
	long double (*survivalProb)(long double, double *), 
	double *survivalProbParameters,
	double Eres)
{
	double binsize = (Ebinmax - Ebinmin) / numberOfBins;
	/*Numeric precission*/
	int N = 10000;

	/*Analysis energy range of THE TRUE electron recoil*/
	long double Tanalmin = 0.001;
	long double Tanalmax = 20;
	long double h = (Tanalmax - Tanalmin) / N;

	long double *countRate = new long double [N+1];

	for(int i = 0; i <= N; i++){
		printf("Computando, %f%...\n", 100*(float)(i+1)/(N+1));
		countRate[i] = countRate_nu_e_Scattering(Flux, Xsection, XsectionParameters, survivalProb, survivalProbParameters, Tanalmin + i*h);
	}

	long double E = Ebinmin;
	for(int i = 0; i < numberOfBins; i++){
		Bins[i] = 0;
		Bins[i] += countRate[0]*resolutionFactor(E, E + binsize, Tanalmin, Eres);
		Bins[i] += countRate[N]*resolutionFactor(E, E + binsize, Tanalmin + N*h, Eres);
		for(int j = 0; j < N; j+=2){
			Bins[i] += 4*countRate[j]*resolutionFactor(E, E + binsize, Tanalmin + j*h, Eres);
			Bins[i] += 2*countRate[j+1]*resolutionFactor(E, E + binsize, Tanalmin + (j+1)*h, Eres);
		}
		Bins[i] *= h/3;
		Bins[i] *= Ntau(E);
		E += binsize;
	}
	delete [] countRate;
}