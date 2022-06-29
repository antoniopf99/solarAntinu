#include"constant_matter_evolution.h"

void pauli_vector::seta_dot_sigma(){
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 2; j++)
			dot_sigma[i][j] = 0;

	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 2; j++)
			for(int k = 0; k < 2; k++)
				dot_sigma[j][k] += n[i]*PAULI[i][j][k];
}

nu_state constant_matter_evolution_2gen(double E, double L, double Ne, double theta, double Dm, nu_state nu_ini){
	double s2 = sin(2*theta);
	double c2 = cos(2*theta);

	double Acc = ((7.6324662178684095e-8)*Ne*2*E); //E is in MeV!
	double Dm_matter = sqrt(pow(Dm*c2 - Acc, 2) + pow(Dm*s2, 2));
	double c2M = (Dm*c2 - Acc)/Dm_matter;
	double s2M = Dm*s2/Dm_matter;

	double phix = 1.26693*(Dm_matter*L)/E;

	pauli_vector pauli;
	pauli.n[0] = s2M;
	pauli.n[1] = 0;
	pauli.n[2] = -c2M;
	pauli.seta_dot_sigma();

	complex<double> evolution[2][2];
	for(int i = 0; i < 2; i++){
		for(int j = 0; j < 2; j++){
			if(i == j)
				evolution[i][i] = cos(phix);
			else
				evolution[i][j] = 0;
		}
	}
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 2; j++)
			evolution[i][j] += -I*sin(phix)*pauli.dot_sigma[i][j];

	nu_state nu_final;
	for(int i = 0; i < 2; i++)
		nu_final.index[i] = 0;

	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 2; j++)
			nu_final.index[i] += evolution[i][j]*nu_ini.index[j];

	return(nu_final);
}

double constant_matter_2gen(double E, double L, double theta, double Dm, double ne){
	nu_state nu;
	nu.index[0] = 1;
	nu.index[1] = 0;

	nu = constant_matter_evolution_2gen(E, L, ne, theta, Dm, nu);

	return(norm(nu.index[0]));
}

double castle_wall_prob_2gen(double E,
	double L,
	double theta,
	double Dm,
	double T1,
	double T2,
	double n1,
	double n2)
{
	nu_state nu_atual;
	nu_atual.index[0] = 1;
	nu_atual.index[1] = 0;

	int periodos = floor(L/T2);

	for(int i = 0; i < periodos; i++){ /*evolving over periods*/
		nu_atual = constant_matter_evolution_2gen(E, T1, n1, theta, Dm, nu_atual); /*first density*/
		nu_atual = constant_matter_evolution_2gen(E, T2-T1, n2, theta, Dm, nu_atual); /*second density*/
	}

	if(L-periodos*T2 > T1){
		nu_atual = constant_matter_evolution_2gen(E, T1, n1, theta, Dm, nu_atual);
		nu_atual = constant_matter_evolution_2gen(E, L-periodos*T2-T1, n2, theta, Dm, nu_atual);
	}
	else{
		nu_atual = constant_matter_evolution_2gen(E, L-periodos*T2, n1, theta, Dm, nu_atual);
	}

	return(norm(nu_atual.index[0]));
}

nu_state constant_matter_evolution_3gen(double E, double L, double Ne, nu_state nu_ini, osc_param op){
	/*the vacuum portion H_0*/
	double H_0[3][3];
	double cte;
	cte = 2533.8653588396974; /* 1e-12 * 1/(2*hc) in units of 1/MeVkm */
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			H_0[i][j] = 0.0;
	H_0[1][1] = cte*L*(op.Dmq21)/E;
	H_0[2][2] = cte*L*(op.Dmq31)/E;

	/*the matter potential, in units of sqrt(2)*G_F*n_e*/
	complex<double> V_matter[3][3];
	double cte2, V_flavor[3][3];
	cte2 = 0.00038679283503942005; /* 1e-6 sqrt(2)*GF */
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			V_matter[i][j] = 0.0;
			V_flavor[i][j] = cte2*L*Ne*op.V[i][j];
		}
	}

	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			V_matter[i][j] += conj(op.U[0][i])*V_flavor[0][0]*(op.U[0][j]);

	/*the full Hamiltonian*/
	complex<double> H[3][3];
	for(int i = 0; i < 3; i++){
		for(int j = 0; j < 3; j++){
			H[i][j] = -(H_0[i][j] + V_matter[i][j]);
			//if(isnan(real(H[i][j])) || isnan(imag(H[i][j]))){
			//	cout << "HAMILTONIANO NAN!" << endl;
			//	cout << "L = " << L << endl;
			//	cout << "Ne = " << Ne << endl;
			//	cout << "E = " << E << endl;
			//}
		}
	}

	//double checa_unitario;
	//for(int i = 0; i < 3; i++){
	//	for(int j = 0; j < 3; j++){
	//		checa_unitario = norm(H[i][j] - conj(H[j][i]));
	//		if(checa_unitario > 0.000001)
	//			cout << "HAMILTONIANO NÃO UNITARIO NA " << i << " " << j << " " << "POR " << checa_unitario << endl;
	//	}
	//}

	nu_state nu_final;

	/*exponential with Putzer's algorithm*/	
	matrixexp_times_vector_3x3_putzer(H, nu_ini.index, nu_final.index);

	//double checa_normalizado_ini = 0, checa_normalizado_fim = 0;
	//for(int i = 0; i < 3; i++){
	//	checa_normalizado_ini += norm(nu_ini.index[i]);
	//	checa_normalizado_fim += norm(nu_final.index[i]);
	//}
	//if(abs(checa_normalizado_fim - checa_normalizado_ini) > 0.1){
	//	cout << "A DIFERENÇA DAS NORMAS DEU " << abs(checa_normalizado_fim - checa_normalizado_ini) << endl;
	//	for(int i = 0; i < 3; i++){
	//		for(int j = 0; j < 3; j++)
	//			cout << H[i][j] << " ";
	//		cout << endl;
	//	}
	//	cout << endl;
	//}

	return(nu_final);
}

double constant_matter_prob_3gen(int nu_emitted, int nu_detected, double E, double L, double Ne, osc_param op){
	/*Initial condition*/
	nu_state nu;
	for(int i = 0; i < 3; i++)
		nu.index[i] = op.U[nu_emitted][i];

	nu = constant_matter_evolution_3gen(E, L, Ne, nu, op);

	/*flavor basis amplitude*/
	complex<double> amp = 0.0;
	for(int i = 0; i < 3; i++)
		amp += conj(op.U[nu_detected][i])*nu.index[i];

	return(norm(amp));
}

double survival_prob_sun_ad_earth_generic(double Rfrac, double E, osc_param op, nu_state (*evolution)(double, nu_state, osc_param)){
	nu_state nu_ini, nu_decoh, nu_earth_evolved;
	complex<double> amp_e;
	double Pe[3] = {0, 0, 0}, eigenv[3];

	/*adiabatic decoherent evolution inside the sun*/
	nu_decoh = sun_adiabatic_decoherent_evolution(Rfrac, E, op, eigenv);


	/*now we evolve the 1 and 2 matter states inside the earth (and use unitarity for the third)*/
	for(int i = 0; i < 2; i++){

		/*initial condition for earth evolution, matter eigenstate*/
		for(int j = 0; j < 3; j++){
			if(j == i)
				nu_ini.index[j] = 1.0;
			else
				nu_ini.index[j] = 0.0;
		}

		/*earth evolution*/
		nu_earth_evolved = evolution(E, nu_ini, op);
			
		/*survival probability*/
		amp_e = 0.0;
		for(int j = 0; j < 3; j++)
			amp_e += op.U[0][j]*nu_earth_evolved.index[j];
		Pe[i] = norm(amp_e);
	}

	Pe[2] = (1 - Pe[0] - Pe[1]);

	/*the survival probs with relative admixtures*/
	return(norm(nu_decoh.index[0])*Pe[0] + norm(nu_decoh.index[1])*Pe[1] + norm(nu_decoh.index[2])*Pe[2]);
}

nu_state earth_evolution_nlayers(double E, nu_state nu_ini, osc_param op, int n_layers, double radius[], double Ne[]){
	nu_state nu_final;

	double L, thetaz = acos(op.cosz);
	int final_layer = 0;

	if(thetaz >= tangent_layer_angle(radius[1], op.h)){
		L = d_outer_single(op.cosz, op.h);
		if(isnan(L))
			cout << "Tá nan na camada " << 0 << endl;
		nu_final = constant_matter_evolution_3gen(E, L, Ne[0], nu_ini, op);
	}
	else{
		L = d_outer_far(op.cosz, radius[1], op.h);
		nu_final = constant_matter_evolution_3gen(E, L, Ne[0], nu_ini, op);
		
		for(int i = 2; i < n_layers && !final_layer; i++){
			if(thetaz < tangent_layer_angle(radius[i], op.h)){
				L = d_inner_2tracks(op.cosz, radius[i-1], radius[i], op.h);
				nu_final = constant_matter_evolution_3gen(E, L, Ne[i-1], nu_final, op);
			}
			else{
				final_layer = i-1;
			}
		}
		if(!final_layer)
			final_layer = n_layers-1;
		
		L = d_inner_single(op.cosz, radius[final_layer], op.h);
		nu_final = constant_matter_evolution_3gen(E, L, Ne[final_layer], nu_final, op);
		
		for(int i = final_layer-1; i >= 1; i--){
			L = d_inner_2tracks(op.cosz, radius[i], radius[i+1], op.h);
			nu_final = constant_matter_evolution_3gen(E, L, Ne[i], nu_final, op);
		}
		
		L = d_outer_near(op.cosz, radius[1], op.h);
		nu_final = constant_matter_evolution_3gen(E, L, Ne[0], nu_final, op);
	}
	return(nu_final);
}

nu_state earth_evolution_2layers(double E, nu_state nu_ini, osc_param op){
	double radius[2], Ne[2];

	radius[0] = R_EARTH, Ne[0] = 0.495*4.5; //Mantle
	radius[1] = 3480, Ne[1] = 0.475*11.5; //Core

	nu_state nu_final = earth_evolution_nlayers(E, nu_ini, op, 2, radius, Ne);

	return(nu_final);
}

double sun_earth_prob_2layers(double Rfrac, double E, osc_param op){	
	if(op.cosz <= 0)
		return(survival_prob_sun_general_adiabatic(Rfrac, E, op));
	return(survival_prob_sun_ad_earth_generic(Rfrac, E, op, &earth_evolution_2layers));
}

double sun_earth_2layers_Average8B(double E, osc_param op){
	/*points used to compose the average of boron neutrino production points*/
	int N = 10;
	double resultado = 0, normalizacao = 0;
	double RfracMax = 0.15;
	double h = RfracMax/N;
	double boronFrac;
	for(double Rfrac = 0; Rfrac < RfracMax; Rfrac += h){
		boronFrac = interpola_dados_solares(Rfrac, 8);
		resultado += boronFrac*sun_earth_prob_2layers(Rfrac, E, op);
		normalizacao += boronFrac;
	}
	return(resultado/normalizacao);
}

nu_state earth_evolution_5layers(double E, nu_state nu_ini, osc_param op){
	double radius[5], Ne[5];

	radius[0] = R_EARTH, Ne[0] = 1.69; //Crust
	radius[1] = 0.937*R_EARTH, Ne[1] = 1.92; //Upper mantle
	radius[2] = 0.895*R_EARTH, Ne[2] = 2.47; //Lower mantle
	radius[3] = 0.546*R_EARTH, Ne[3] = 5.24; //Outer core
	radius[4] = 0.192*R_EARTH, Ne[4] = 6.05; //Inner core

	nu_state nu_final = earth_evolution_nlayers(E, nu_ini, op, 5, radius, Ne);

	return(nu_final);
}

double sun_earth_prob_5layers(double Rfrac, double E, osc_param op){
	if(op.cosz <= 0)
		return(survival_prob_sun_general_adiabatic(Rfrac, E, op));
	return(survival_prob_sun_ad_earth_generic(Rfrac, E, op, &earth_evolution_5layers));
}

double atmospheric_neutrino_5layers(int flavor_ini, int flavor_end, double E, osc_param op){
	nu_state nu_ini, nu_final;
	
	for(int i = 0; i < 3; i++)
		nu_ini.index[i] = op.U[flavor_ini][i]; //Inital condition

	double a = 15;

	double cosz = op.cosz;
	double tan2 = (1 - cosz*cosz)/(cosz*cosz);
	double L_atm = a*abs(cosz)*(1 + (2*tan2)/(sqrt(1 + 2*(a + op.h)*tan2/R_EARTH) + sqrt(1 + 2*(op.h)*tan2/R_EARTH)));
	nu_final = constant_matter_evolution_3gen(E, L_atm, 0, nu_ini, op);
	
	nu_final = earth_evolution_5layers(E, nu_final, op);

	complex<double> amp = 0;

	for(int i = 0; i < 3; i++)
		amp += op.U[flavor_end][i]*nu_final.index[i];

	return(norm(amp));
}