#ifndef SOLARFLUX_H
#define SOLARFLUX_H

/*Returns solar continuous neutrino flux in cm⁻²s⁻¹keV⁻¹
type = 0 -> all continuous neutrino fluxes
type = 1 -> pp flux
type = 2 -> hep flux
type = 3 -> 8B flux
type = 4 -> 13N flux
type = 5 -> 15O flux
type = 6 -> 17F flux*/
long double solarFlux(int type, long double Enu);

#endif