#ifndef EARTH_MODEL_H
#define EARTH_MODEL_H

/*
Earth's radius in km
*/
const double R_EARTH = 6371.0; 

/*
Returns the density of the planet as a function of the fractional radius,
in units of g/cm³
*/
double earth_density(double x/*Rfrac*/);
double earth_density_5layers(double x/*Rfrac*/);

/*
Returns the electron density of the planet as a function of the fractional radius,
in units of NA/cm³
*/
double earth_electron_density_NA(double Rfrac);

/*
Returns the distance a neutrino travels inside the planet as a function of the cossine
of the solar zenith angle, in units of R_EARTH
*/
double distance_inside_earth(double cosz, double h);

/*
Returns the distance a neutrino travels inside the an internal layer of the planet of
radius r_core (km). Value is returned in km
*/
double distance_inside_core(double r_core, double cosz);

/*
Returns the radial distance from the center of the planet, for a neutrino that already 
traveled dprime along a tragectory given by the cossine of the solar zenith angle cosz,
in units of R_EARTH 
*/
double radial_distance_core(double cosz, double h, double dprime);

/*
Returns the angle that corresponds to a neutrino path that is tangent to a given layer.
Takes into account the detector depth h.
*/
double tangent_layer_angle(double r_layer, double h);

/*
Returns the distance traveled inside each of the earths layers in a generic
n layer model. Takes into account the detector's depth h. See arXiv:2005.07719v2
returns result in units of km!
*/
double d_outer_single(double cosz, double h);
double d_outer_near(double cosz, double r, double h);
double d_outer_far(double cosz, double r, double h);
double d_inner_single(double cosz, double r, double h);
double d_inner_2tracks(double cosz, double r_bigger, double r_smaller, double h);


#endif