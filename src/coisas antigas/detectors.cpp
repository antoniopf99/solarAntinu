#include"detectors.h"

long double exposureXsection(long double E){
	double year = 3.1536; /*year in seconds, FACTOR OF e7*/
	double pot = 4.56; /*protons on target, FACTOR OF e31*/
	/*the cross section has a factor of e-42 missing, so we need e-42+7+31 = e-4*/
	return(year * pot * pow(10,-4) * inverseBeta(E));
}

double detectorEfficiency(long double x, const string datafile){
	double x1, x2, y1, y2, resultado;
	bool foi = false;

	/*********************/
	/*DATA DIRECTORY NAME*/
	/*********************/
	string datadir = "data/";

	if(datafile == "")
		resultado = 85;
	else{
		ifstream arquivo;
		arquivo.open(datadir + datafile);
		arquivo >> x1 >> y1;
		if(x1 > x){
 			resultado = y1;
 			foi = true;
 		}
		while(!arquivo.eof() && !foi){
 			arquivo >> x2 >> y2;
 			if(x1 <= x && x <= x2){
 				resultado = ((y2-y1)/(x2-x1))*(x-x1) + y1;
 			 	foi = true;
 			}
 			x1 = x2;
 			y1 = y2;
 			arquivo.peek();
 		}
 		arquivo.close();
 		if(!foi)
 			resultado = y1;
 	}
 	return(resultado/100);
}

long double energyResolution(double Eobs, double Ereal, double Res){
	if(Ereal <= 0){
		cout << "ERROR COMPUTING EVENT AT ZERO ENERGY! (At the energyResolution function)." << endl;
		return(0);
	}
	long double d = ( (Res/100) * sqrt(Ereal) ) / (2 * sqrt(2*log(2)));
	long double N = 1 / (sqrt(2*M_PI) * d);
	return(N * exp(-0.5 * pow((Eobs - Ereal)/d,2)));
}

long double resolutionFactor(double Eini, double Efinal, double Ee, double Res, int numericIntervals){
	long N = 2 * numericIntervals;/*number of sub divisions must be even!*/
	long double h = (Efinal - Eini) / N;
	long double soma = energyResolution(Eini,Ee,Res) + energyResolution(Eini + N*h,Ee,Res);
	for(long i = 1; i < N; i+=2)
		soma += 4*energyResolution(Eini + i*h,Ee,Res) + 2*energyResolution(Eini + (i+1)*h,Ee,Res);
	soma *= h/3;
	return(soma);
}

void zeraBins(double Ntheo[][3], int numberOfBins){
	for(int i = 0; i < numberOfBins; i++)
		Ntheo[i][0] = 0;
}

void calculaBins(vector<Reactor> &reatores, long double distance, double res, const string datafile, long numericIntervals, double Ntheo[][3], double BinSize, int numberOfBins, double EanalMin, double EanalMax, bool oscila, long double theta12, long double theta13, long double Dmq){
 	double Bini = EanalMin, Bfim;
 	vector<long double> Fluxos_etc;
 	Fluxos_etc.reserve(numericIntervals+2);
 	double fluxoparticular;
 	long double h = (EanalMax - EanalMin) / numericIntervals;
 	for(long i = 0; i <= numericIntervals; i++){
 		fluxoparticular = 0;
 		for(Reactor it : reatores){
 			if(it.distance < distance){
 				if(oscila){
 					fluxoparticular += it.Fluxo(EanalMin + i*h) * survivalProbConstMatter(it.distance * pow(10,-5)/*cm to km*/,(EanalMin + i*h) * pow(10,-3)/*MeV to GeV*/, Dmq, theta12, theta13);
 				}
 				else
 					fluxoparticular += it.Fluxo(EanalMin + i*h);
 			}
 		}
 		fluxoparticular *= exposureXsection(EanalMin + i*h) * detectorEfficiency(EanalMin + i*h - 0.7823, datafile); 
 		Fluxos_etc[i] = fluxoparticular;
 	}
	long double pinf, psup, adiciona;
 	for(int i = 0; i < numberOfBins; i++){
 		Bfim = Bini + BinSize;
 		adiciona = Fluxos_etc[0] * resolutionFactor(Bini - 0.7823, Bfim - 0.7823, EanalMin - 0.7823, res, 2*(numericIntervals/40)) + Fluxos_etc[numericIntervals] * resolutionFactor(Bini - 0.7823, Bfim - 0.7823, EanalMin + numericIntervals*h - 0.7823, res, 2*(numericIntervals/40));
 		for(long j = 1; j < numericIntervals; j+=2){
			pinf = Fluxos_etc[j] * resolutionFactor(Bini - 0.7823, Bfim - 0.7823, EanalMin + j*h - 0.7823, res, 2*(numericIntervals/40));
			psup = Fluxos_etc[j+1] * resolutionFactor(Bini - 0.7823, Bfim - 0.7823, EanalMin + (j + 1)*h - 0.7823, res, 2*(numericIntervals/40));
			adiciona += 4*pinf + 2*psup;
 		}
 		adiciona *= h/3;
 		Ntheo[i][0] += adiciona;
 		Ntheo[i][1] = Bini;
 		Ntheo[i][2] = Bfim;
		/*cout << "Bin number " << i + 1 << " done." << endl;*/
 		Bini = Bfim;
 	}
}