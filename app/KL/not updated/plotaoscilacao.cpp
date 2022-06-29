#include<iostream>
#include<fstream>
#include<cstring>
#include"oscilation.h"

using namespace std;

int main(){

	float PASSO = 0.0001;
	
	string outdir = "results/"; 
	
	int mais;
	do{
		int inicial, final;
		ofstream saida;
		string nomesaida;
		cout << "Output file name: ";
		cin >> nomesaida;
		saida.open(outdir + nomesaida);
		cout << "Initial neutrino and final neutrino type (0=e, 1=m, 2=t): ";
		cin >> inicial >> final;
		/*10^3.257 keV ~ 1.8 MeV inverse beta decay treshold*/
		for(float pot = 3.257; pot < 4; pot += PASSO)
			saida << pow(10,pot) << " " << survivalProb(inicial,final,180,pow(10,-6)*pow(10,pot),true) << endl;
		/*
		for(int E = pow(10,3.257); E < pow(10,4); E += PASSO)
			saida << E << " " << survivalProb(inicial,final,180,pow(10,-6)*E,true) << endl;
		*/
		saida.close();
		cout << "More? (1=y, 0=n): ";
		cin >> mais;
	}while(mais);
	return(0);
}
