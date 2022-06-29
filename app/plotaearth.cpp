#include<iostream>
#include<fstream>
#include"crosssections.h"
#include"earth_model.h"

using namespace std;

int main(){
	ofstream saida;
	saida.open("results/earthdensitypath.txt");

	zprime_param zp;
	zp.antinu = false;
	zp.gX = zp.a = zp.e = zp.mZprime = zp.QlL = zp.QlR = zp.QnuL = 0;

	zp.Enu = 10;

	for(zp.T = 0; zp.T < 10; zp.T += 0.01){
		cout << zp.T << " ";
		zp.NC_and_CC = true;
		cout << nueScatteringJUNO(zp) << " ";
		zp.NC_and_CC = false;
		cout << nueScatteringJUNO(zp) << " ";
		cout << 6*nueScatteringJUNO(zp) << endl;
	}


	//for(double Rfrac = 0; Rfrac < 1; Rfrac += 0.001){
	//	cout << Rfrac << " ";
	//	cout << earth_electron_density_NA(Rfrac) << " ";
	//	cout << earth_density_5layers(Rfrac) << endl;
	//}

	saida.close();
	return(0);
}