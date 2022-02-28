#include<iostream>
#include<fstream>
#include<cstring>
#include"solarmodel.h"

using namespace std;

int main(){
	ofstream saida;
	string filename, outdir = "results/";
	int more;

	do{
		cout << "Enter with output file name: ";
		cin >> filename;

		saida.open(outdir + filename);

		seta_dados_Ne();

		//for(double Rfrac = 0; Rfrac < 1; Rfrac += 0.001){
		//	saida << Rfrac << " " << electron_density_sun_NA(Rfrac) << endl;
		//	//saida << Rfrac << " " << interpola_electron_density(Rfrac) << endl;
		//}
		
		int dado;
		cout << "Enter with the data you want: ";
		cin >> dado;
		double h = 0.001, N = 1;
		//For normalizing the fraction of neutrino production
		N = 0;
		for(long double Rfrac = 0; Rfrac < 0.5; Rfrac += h){
			N += h * interpola_dados_solares(Rfrac, dado);
		}
		for(long double Rfrac = 0; Rfrac < 0.5; Rfrac += h){
			saida << Rfrac << " " << interpola_dados_solares(Rfrac, dado) / N << endl;
		}

		saida.close();
		cout << "Do you want more data? ";
		cin >> more;
	}while(more);
	return(0);
}