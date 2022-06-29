#include"crosssections.h"

double atanh(double x){
	return(0.5*(log(x+1) - log(x-1)));
}

double inverseBeta(double E){
	if(E < 1.805)
		return(0);
	double Me = 0.511;
	double Ee = E - 1.293;
	return(pow(10,-43+42) * sqrt(pow(Ee,2) - pow(Me,2)) * Ee * pow(E,-0.07056 + 0.02018*log(E)-0.001953*pow(log(E),3)));
	//return(9.75 * 3.154 * pow(10,-3/*-42+7+32*/) * pow(E/10,2));
}

double nueScatteringSM(double Enu, double T){
	/**************************/
	/*CROSS SECTION PARAMETERS*/
	/**************************/
	double GF = 1.1663787e-11; /*Fermi constant in MeV-2*/
	double me = 0.5109989461; /*Electron mass in MeV*/
	double hc = 1.973269804e-11; /*hc in units of MeV cm*/
	
	double No = hc*hc*GF*GF*me/(4*M_PI);

	double sw2 = 0.22290; /*Sine of Weinberg angle squared*/
	double g_2 = -2*sqrt(2)*sw2;
	double g_1 = -sqrt(2)*(1+2*sw2);

	double Tmax = 2*pow(Enu,2) / (me + 2*Enu);
	double Enumin = 0.5*(T+sqrt(T*(T+2*me)));

	if(T > Tmax || Enu < Enumin)
		return(0);

	return(No * (pow(g_2,2) + pow(g_1,2)*pow(1 - T/Enu,2) - g_1*g_2*((me*T)/pow(Enu,2))));
}

double nueScatteringJUNO(zprime_param p){
	double s0 = 88.06e-46;
	double me = 0.511;
	double g1, g2;
	if(p.NC_and_CC){
		g1 = 0.73;
		g2 = 0.23;
	}
	else{
		g1 = -0.27;
		g2 = 0.23;
	}
	return((s0/me)*(g1*g1 + g2*g2*pow(1 - p.T/p.Enu, 2) - g1*g2*me*p.T/pow(p.Enu, 2)));
}

double nueScatteringBSM(zprime_param p){
	/**************************/
	/*CROSS SECTION PARAMETERS*/
	/**************************/
	
	/*fixed*/
	double GF = 1.1663787e-11; /*Fermi constant in MeV-2*/
	double me = 0.5109989461; /*Electron mass in MeV*/
	double hc = 1.973269804e-11; /*hc in units of MeV cm*/
	
	double No = hc*hc*GF*GF*me/(4*M_PI);

	double sw2 = 0.22290; /*sine of Weinberg angle squared*/
	double sw = sqrt(sw2);
	double cw2 = 1 - sw2;
	double cw = sqrt(cw2);
	double g = 0.652905; /*coupling to the W*/
	double g2 = g*g;

	double Enu = p.Enu;
	double T = p.T;

	/*model dependent*/
	double mZprime = p.mZprime;

	double ca = cos(p.a);
	double sa = sin(p.a);
	double ca2 = ca*ca;
	double sa2 = sa*sa;

	double gX = p.gX;
	double gX2 = gX*gX;
	double e = p.e;
	double e2 = e*e;
	double sr2 = sqrt(2);

	double g1SM = -2*sr2*sw2*ca2;
	double g2SM = sr2*ca2*(1-2*sw2)-(p.NC_and_CC)*2*sr2;

	double QlL = p.QlL; 
	double QlR = p.QlR;
	double QnuL = p.QnuL; 

	double a1 = sa2*( gX2*((sr2*cw2*QnuL*QlR)/(g2*(e2-1))) + e*gX*((sr2*sw*cw*(2*QnuL+QlR))/(g*(e2-1))) + e2*((2*sr2*sw2)/(e2-1)) ) + ca*sa*( -gX*((sr2*sqrt(1-e2)*cw*(2*sw2*QnuL  +  QlR))/(g*(e2-1))) - e*((2*sr2*sqrt(1-e2)*sw*(sw2+1))/(e2-1)) );
	double a2 = sa2*( gX2*((sr2*cw2*QnuL*QlL)/(g2*(e2-1))) + e*gX*((sr2*sw*cw*(QlL + QnuL))/(g*(e2-1))) + e2*((2*sr2*sw2)/(e2-1)) ) + ca*sa*( -gX*((sr2*sqrt(1-e2)*cw*((2*sw2-1)*QnuL+QlL))/(g*(e2-1))) - e*((2*sr2*sqrt(1-e2)*sw2* sw)/(e2-1)) );
	double b1 = ca2*( (g2*e2*sw2)/(2*(e2-1)*cw2) + (g*e*gX*sw*QnuL)/(2*(e2-1)*cw) + (g*e*gX*sw*QlR)/(4*(e2-1)*cw) + (gX2*QnuL*QlR)/(4*(e2-1)) ) + ca*sa*( e*((g2*sqrt(1-e2)*sw*(sw2+1))/(2*(e2-1)*cw2)) + gX*((g*sqrt(1-e2)*(2*sw2*QnuL+QlR))/(4*(e2-1)*cw)) ) + sa2*(-(g2*sw2)/(2*cw2)); 
	double b2 = ca2*( (g2*e2*sw2)/(4*(e2-1)*cw2) + (g*e*gX*sw*(QlL+QnuL))/(4*(e2-1)*cw) + (gX2*QnuL*QlL)/(4*(e2-1)) ) + ca*sa*( -e*((g2*sw2*sw)/(2*sqrt(1-e2)*cw2)) - gX*( (g*((2*sw2-1)*QnuL + QlL) ) /(4*sqrt(1-e2)*cw)) ) + sa2*((g2*(1-2*sw2))/(4*cw2));

	double r = (1 / ((2*me*T+pow(mZprime,2))*GF));

	double Tmax = 2*pow(Enu,2) / (me + 2*Enu);
	double Enumin = 0.5*(T + sqrt(T*(T+2*me)));

	if(T > Tmax || Enu < Enumin)
		return(0);

	double g_1;
	double g_2;
	if(p.antinu){
		g_1 = g1SM + a1 + b1*r;
		g_2 = g2SM + a2 + b2*r;
	}
	else{
		g_1 = g2SM + a2 + b2*r;
		g_2 = g1SM + a1 + b1*r;
	}
	return(No * (pow(g_1,2) + pow(g_2,2)*pow(1 - T/Enu,2) - g_1*g_2*((me*T)/pow(Enu,2))));
}

double partiaNueXsection(bool antinu, double gX, double a, double e, double mZprime, double QlL, double QlR, double QnuL, double Enu, double T1, double T2){
	/**************************/
	/*CROSS SECTION PARAMETERS*/
	/**************************/
	
	/*fixed*/
	double GF = 1.166378 * pow(10,-11); /*Fermi constant in MeV^-2*/
	double me = 0.510998; /*electron mass in MeV*/
	double sw2 = 0.22290; /*sine of Weinberg angle squared*/
	double sw = sqrt(sw2);
	double cw2 = 1 - sw2;
	double cw = sqrt(cw2);
	double g = 0.652905; /*coupling to the W*/

	double gX2 = pow(gX,2);
	double e2 = pow(e,2);
	double g2 = pow(g,2);
	double sr2 = sqrt(2);

	double ca = cos(a);
	double sa = sin(a);
	double ca2 = pow(ca,2);
	double sa2 = pow(sa,2);
	
	/*model dependent*/
	double g1SM = -2*sr2*sw2*ca2;
	double g2SM = sr2*ca2*(1-2*sw2)-2*sr2;

	double Aplus = (g*e*cw*gX*sw*(2*QnuL + QlR) + cw2*gX2*QnuL*QlR + 2*pow(g*e,2)*sw2) / (2*(1 - e2));
	double Bplus = (-g*sqrt(1 - e2)*cw*gX*(2*sw2*QnuL + QlR) - 2*g2*e*sqrt(1 - e2)*pow(sw,3) - 2*g2)*e*sqrt(1 - e2)*sw / (2*(1 - e2));
	double Aminus = ((-cw*gX*QlL - g*e*sw) * (cw*gX*QnuL + g*e*sw)) / (2*(1 - e2));
	double Bminus = (g*sqrt(1 - e2)*(cw*gX*(-QnuL + 2*sw2*QnuL + QlL) + 2*g*e*pow(sw,3))) / (2*(1 - e2));

	double a1 = (-2*sr2*(Aplus*sa2 + Bplus*ca*sa)) / g2;
	double b1 = -(g2*sa2*sw2 + Aplus*ca2 - Bplus*ca*sa) / (2*cw2);
	double a2 = (-2*sr2*(Aminus*sa2 + Bminus*ca*sa)) / g2;
	double b2 = -(g2*sa2*(sw2 - 0.5) + Aminus*ca2 - Bminus*ca*sa) / (2*cw2);

	double x[4];

	if(antinu){
		x[0] = g2SM + a2;
		x[1] = g1SM + a1;
		x[2] = b2;
		x[3] = b1;
	}
	else{
		x[0] = g1SM + a1;
		x[1] = g2SM + a2;
		x[2] = b1;
		x[3] = b2;
	}

	double I[4][4];

	for(int i = 0; i < 4; i++)
		for(int j = 0; j < 4; j++)
			I[i][j] = 0;

	I[0][0] = -(me*(T1-T2)*(3*pow(Enu,2) - 3*(T1+T2)*Enu + pow(T1,2) + pow(T2,2) + T1*T2))/(12*M_PI*pow(Enu,2));	
	I[0][1] = (pow(me,2)*(T1-T2)*(T1+T2)) / (8*M_PI*pow(Enu,2));
	I[0][2] = (log((pow(mZprime,2) + 2*me*T2) / (pow(mZprime,2) + 2*me*T1))*pow(pow(mZprime,2) + 2*Enu*me,2) + 2*me*(pow(mZprime,2) + me*(4*Enu - T1 - T2))*(T1 - T2)) / (16*M_PI*pow(Enu,2)*GF*pow(me,2));
	I[1][1] = (me*(T2 - T1)) / (4*M_PI);
	I[1][2] = -((atanh((me*(T1 - T2)) / (pow(mZprime,2) + me*(T1 + T2)))*pow(mZprime,2) + me*(T2 - T1)) / (8*M_PI*pow(Enu,2)*GF));
	I[1][3] = log((pow(mZprime,2) + 2*me*T2) / (pow(mZprime,2) + 2*me*T1)) / (4*M_PI*GF);
	I[2][2] = -((0.5*log((pow(mZprime,2) + 2*me*T2) / (pow(mZprime,2) + 2*me*T1))*(pow(mZprime,2) + 2*Enu*me) + (me*(T1 - T2)*(pow(mZprime,4) + me*(2*Enu + T1 + T2)*pow(mZprime,2) + 2*pow(me,2)*(pow(Enu,2) + T1*T2)))/((pow(mZprime,2) + 2*me*T1)*(pow(mZprime,2) + 2*me*T2))) / (8*M_PI*pow(Enu,2)*pow(GF,2)*pow(me,2)));
	I[2][3] = -(((1/(pow(mZprime,2) + 2*me*T2) - 1/(pow(mZprime,2) + 2*me*T1))*pow(mZprime,2) - log(pow(mZprime,2) + 2*me*T1) + log(pow(mZprime,2) + 2*me*T2)) / (16*M_PI*pow(Enu,2)*pow(GF,2)));
	I[3][3] = (me*(T2 - T1)) / (4*M_PI*pow(GF,2)*(pow(mZprime,2) + 2*me*T1)*(pow(mZprime,2) + 2*me*T2));

	double result = 0;

	for(int i = 0; i < 4; i++)
		for(int j = 0; j < 4; j++)
			result += x[i]*I[i][j]*x[j];
	result *= pow(GF,2);

	return(result);
}