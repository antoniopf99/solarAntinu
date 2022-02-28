#include"geoCoords.h"

long double Distancia(geoCoord loc1, geoCoord loc2){

	/**********************/
	/*EARTH'S RADIUS IN CM*/
	/**********************/
	long double RT = 6.371 * pow(10,8);

	double theta1, phi1, theta2, phi2, angulo;
	theta1 = emRadianos(loc1.lat + 90);
	phi1 = emRadianos(loc1.lon);
	theta2 = emRadianos(loc2.lat + 90);
	phi2 = emRadianos(loc2.lon);
	angulo = acos(sin(theta1)*sin(theta2)*cos(phi1)*cos(phi2)+sin(theta1)*sin(theta2)*sin(phi1)*sin(phi2)+cos(theta1)*cos(theta2));
	return(2 * RT * sin(angulo / 2));
}

long double emRadianos(long double graus){
	return(graus / 57.2958);
}
