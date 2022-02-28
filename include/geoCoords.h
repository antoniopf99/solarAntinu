#ifndef GEOCOORDS_H
#define GEOCOORDS_H
#include<cmath>

/*COORDINATES*/
typedef struct coord{
	long double lat, lon;
} geoCoord;

/*CONVERT TO RADIANS*/
long double emRadianos(long double graus);

/*COMPUTES DISTANCE BETWEEN TWO COORDINATES(SPHERICAL EARTH)*/
long double Distancia(geoCoord loc1, geoCoord loc2);

#endif