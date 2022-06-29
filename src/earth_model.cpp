#include<cmath>
#include"earth_model.h"

/*A. M. Dziewonski and D. L. Anderson. Preliminary reference earth model. Phys. Earth Planet.
Interiors, 25:297–356, 1981*/
double earth_density(double x/*Rfrac of the earth*/){
	double result;
	if(x < (1221.5)/R_EARTH){ /*Inner core*/
		result = 13.0885 - 8.8301*x*x;
	}
	else if(x < (3480.0/R_EARTH)){ /*Outer core*/
		result = 12.5815 - 1.2638*x - 3.6426*x*x - 5.5281*x*x*x;
	}
	else if(x < (5701.0/R_EARTH)){ /*Lower mantle*/
		result = 7.9565 - 6.4761*x + 5.5283*x*x - 3.0807*x*x*x;
	}
	else if(x < (5771.0/R_EARTH)){ /*Transition zone I*/
		result = 5.3197 - 1.4836*x;
	}
	else if(x < (5971.0/R_EARTH)){ /*Transition zone II*/
		result = 11.2494 - 8.0298*x;
	}
	else if(x < (6151.0/R_EARTH)){ /*Transition zone III*/
		result = 7.1089 - 3.8045*x;
	}
	else if(x < (6346.6/R_EARTH)){ /*LVZ and LID*/
		result = 2.6910 + 0.6924*x;
	}
	else if(x < (6356.0/R_EARTH)){ /*Crust I*/
		result = 2.900;
	}
	else if(x < (6368.0/R_EARTH)){ /*Crust II*/
		result = 2.600;
	}
	else if(x < 1){ /*Ocean*/
		result = 1.020;
	}
	else{
		result = 0;
	}
	return(result);
}

double earth_density_5layers(double x/*Rfrac of the earth*/){
	if(x > 0.937)
		return(1.69);
	if(x > 0.895)
		return(1.92);
	if(x > 0.546)
		return(2.47);
	if(x > 0.192)
		return(5.24);
	return(6.05);
}

double earth_electron_density_NA(double Rfrac){
	if(Rfrac > 0.546)
		return(0.497*earth_density(Rfrac)); /*giunti*/
	return(0.468*earth_density(Rfrac));
}

double distance_inside_earth(double cosz, double h){
	return(cosz + sqrt(cosz*cosz + 2*h/R_EARTH) - (h/R_EARTH)*(cosz + abs(cosz)));
}

double distance_inside_core(double r_core, double cosz){
	if(R_EARTH*cosz < sqrt(R_EARTH*R_EARTH - r_core*r_core))
		return(0);
	double d = cosz*R_EARTH + sqrt(R_EARTH*R_EARTH*(cosz*cosz - 1) + r_core*r_core);
	return(2*d - 2*R_EARTH*cosz);
}

double radial_distance_core(double cosz, double h, double dprime){
	double c = sqrt(1 - cosz*cosz);// sinz = sqrt(1 - cosz*cosz)
	double d = distance_inside_earth(cosz, h);
	return(sqrt(c*c + pow(d/2 - dprime, 2)));
}

double tangent_layer_angle(double r_layer, double h){
	return(asin(r_layer/(R_EARTH - h)));
}

double d_outer_single(double cosz, double h){
	return(R_EARTH*(cosz + sqrt(cosz*cosz + 2*h/R_EARTH) - (h/R_EARTH)*(cosz + abs(cosz))));
}

double d_outer_near(double cosz, double r, double h){
	double cos2 = cosz*cosz;
	double sin2 = 1 - cos2;
	
	double second_term = sqrt(pow(r/R_EARTH, 2) - (1 - 2*h/R_EARTH)*sin2);

	return(R_EARTH*((1 - h/R_EARTH)*cosz - second_term));
}

double d_outer_far(double cosz, double r, double h){
	double cos2 = cosz*cosz;
	double sin2 = 1 - cos2;
	
	double second_term = sqrt(pow(r/R_EARTH, 2) - (1 - 2*h/R_EARTH)*sin2);

	return(R_EARTH*(sqrt(cos2 + (2*h/R_EARTH)*sin2) - second_term));
}

double d_inner_single(double cosz, double r, double h){
	double sin2 = 1 - cosz*cosz;

	return(2*R_EARTH*sqrt(pow(r/R_EARTH, 2) - (1 - 2*h/R_EARTH)*sin2));
}

double d_inner_2tracks_terms(double cosz, double r, double h){
	double sin2 = 1 - cosz*cosz;

	return(sqrt(pow(r/R_EARTH, 2) - (1 - 2*h/R_EARTH)*sin2));
}

double d_inner_2tracks(double cosz, double r_bigger, double r_smaller, double h){
	return(R_EARTH*(d_inner_2tracks_terms(cosz, r_bigger, h) - d_inner_2tracks_terms(cosz, r_smaller, h)));
}