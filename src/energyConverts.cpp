#include<cmath>

long double emMeVs(long double pot){
	return(pot * 6.2415093433 * pow(10,18));
}

long double emKeV(long double energia){
	return(energia * pow(10,3));
}

long double emMeV(long double energia){
	return(energia * pow(10,-3));
}