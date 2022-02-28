#include<iostream>
#include<cmath>
#include<fstream>
#include"oscilation.h"

using namespace std;

int main(){
	double E = 10;
	for(double x = 0; x < 696342; x += 1)
		cout << x << " " << survivalProb(0, 0, x, E, true) << endl;
 	return(0);
}