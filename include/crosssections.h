#ifndef CROSSSECTIONS_H
#define CROSSSECTIONS_H

#include<cmath>

/*Returns the cross section of the inverse beta decay times 10^42,
as a function of the incoming neutrino energy in MeV, in units of cm^2*/
long double inverseBeta(long double E);


/*Returns the cross section of antinue and e scattering,
as a function of the incoming neutrino energy Enu in MeV and the
electron recoil energy T in MeV, in units of MeV^-1cm^2*/
long double antinueScatteringSM(long double Enu, long double T);
/*Returns the cross section of antinue and e scattering times 10^45,
as a function of the incoming neutrino energy Enu in MeV and the
electron recoil energy T in MeV, in units of MeV^-1cm^2*/
long double antinueScatteringBSM(bool antinu, long double gX, long double a, long double e, long double mZprime, double QlL, double QlR, double QnuL, long double Enu, long double T);

#endif