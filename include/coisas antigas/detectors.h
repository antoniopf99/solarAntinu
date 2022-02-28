#ifndef DETECTORS_H
#define DETECTORS_H

#include<iostream>
#include<cmath>
#include<fstream>
#include<vector>
#include<cstring>
#include<string>
#include"reactors.h"
#include"oscilation.h"
#include"crosssections.h"

using namespace std;

/*CROSS SECTION + EXPOSURE + DETECTOR DETAILS*/
long double exposureXsection(long double E);

/*DETECTOR EFFICIENCY FUNCTION (EVENT ENERGY, NAME OF THE FILE CONTAINING RESOLUTION INFORMATION)
the function will perform a linear interpolation to calculate the efficiency at x*/
double detectorEfficiency(long double x, const string datafile);

/*ENERGY RESOLUTION FUNCTION OF THE DETECTOR
It's a Gaussian centered around Ereal with full-width-at-half-maximum sqrt(E) * Res/100
Res is usually expresed in X%/sqrt(E(MeV))*/
long double energyResolution(double Eobs, double Ereal, double Res);

/*ENERGY RESOLUTION FACTOR
Integrates the energy resolution funtion centered at Ee from Eini to Efinal,
using 1/3 Simpson's with 2*numericInterval sub-divisions*/
long double resolutionFactor(double Eini, double Efinal, double Ee, double Res, int numericIntervals);

/*Insert zeros in every Ntheo[0][] entry*/
void zeraBins(double Ntheo[][3], int numberOfBins);

/*BIN EVENTS CALCULATION
This function adds the events corresponding to the reactors closer than distance in the vector "reactores" in Ntheo[0][i][i+1], of the bin
with lower bound Ntheo[1][i][i+1] and upper bound Ntheo[2][i][i+1], and bin siz BinSize. The numerical integration in
this interval is performed with a 1/3 Simpson's with 2*numericInterval sub-divisions, the integration getting the tail
of events with energies from EanalMin to EanalMax. There should be numberOfBins bins, and the analysis is performed with Dmq and
theta and oscillation parameters. The reactor resolution is res sqrt(MeV)*/
void calculaBins(vector<Reactor> &reatores, long double distance, double res, const string datafile, long numericIntervals, double Ntheo[][3], double BinSize, int numberOfBins, double EanalMin, double EanalMax, bool oscila, long double theta12, long double theta13, long double Dmq);

#endif