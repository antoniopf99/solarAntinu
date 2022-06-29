#include"statistics.h"
#include"detectors.h"

#define MAC p.N_obs[i]&&(p.T[i]||p.BG[i])

double chisquare_twopulls_notminimized(const gsl_vector *x, void *params){
	stat_params p = *(stat_params *)params;
	double ed = gsl_vector_get(x, 0);
	double eb = gsl_vector_get(x, 1);
	double N_pred, result = 0;
	for(int i = 0;  i < p.nBins; i++){
		if(MAC){
			N_pred = (1.0+ed)*p.T[i] + (1.0+eb)*p.BG[i];
			result += 2.0*(N_pred - p.N_obs[i] + p.N_obs[i]*log(p.N_obs[i]/N_pred));
		}
		else{
			N_pred = (1.0+ed)*p.T[i] + (1.0+eb)*p.BG[i];
			result += 2.0*(N_pred - p.N_obs[i] + p.N_obs[i]*1);
		}
	}
	result += pow(ed/p.sd, 2) + pow(eb/p.sb, 2);
	return(result);
}

void chisquare_twopulls_grad(const gsl_vector *x, void *params, gsl_vector *g){
	stat_params p = *(stat_params *)(params);
	double ed = gsl_vector_get(x, 0);
	double eb = gsl_vector_get(x, 1);
	double N_pred, resultd = 0, resultb = 0;

	for(int i = 0; i < p.nBins; i++){
		if(MAC){
			N_pred = (1.0+ed)*p.T[i] + (1.0+eb)*p.BG[i];
			resultd += 2.0*p.T[i]*(1.0 - p.N_obs[i]/N_pred);
			resultb += 2.0*p.BG[i]*(1.0 - p.N_obs[i]/N_pred);
		}
		else{
			N_pred = (1.0+ed)*p.T[i] + (1.0+eb)*p.BG[i];
			//resultd += 2.0*p.T[i]*(1.0 - p.N_obs[i]/N_pred);
			//resultb += 2.0*p.BG[i]*(1.0 - p.N_obs[i]/N_pred);
		}
	}
	resultd += 2.0*ed/pow(p.sd, 2);
	resultb += 2.0*eb/pow(p.sb, 2);

	gsl_vector_set(g, 0, resultd);
	gsl_vector_set(g, 1, resultb);
}

void chisquare_twopulls_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *g){
	*f = chisquare_twopulls_notminimized(x, params);
	chisquare_twopulls_grad(x, params, g);
}

double chisquare_twopulls_JUNO(stat_params p, bool verbo){
	/*providing the function to minimize and its gradient*/
	gsl_multimin_function_fdf chi_2;
	chi_2.n = 2;
	chi_2.f = &chisquare_twopulls_notminimized;
	chi_2.df = &chisquare_twopulls_grad;
	chi_2.fdf = &chisquare_twopulls_fdf;
	chi_2.params = (void *)(&p);

	/*minimization algorithm*/
	const gsl_multimin_fdfminimizer_type *type;
	type = gsl_multimin_fdfminimizer_vector_bfgs2;

	/*starting point*/
	gsl_vector *x = gsl_vector_alloc(2);
	gsl_vector_set(x, 0, 0);
	gsl_vector_set(x, 1, 0);

	gsl_multimin_fdfminimizer *ws = gsl_multimin_fdfminimizer_alloc(type, 2);
	gsl_multimin_fdfminimizer_set(ws, &chi_2, x, 1e-1, 1e-1);

	size_t iter = 0, max_iter = 100;
	int status;

	if(verbo)
		cout << "Começando a minimizar..." << endl;
	do{
		iter++;
		status = gsl_multimin_fdfminimizer_iterate(ws);

		if(status)
			break;

		if(isnan(ws->f))
			break;

		status = gsl_multimin_test_gradient(ws->gradient, 1e-3);
		x = gsl_multimin_fdfminimizer_gradient(ws);

		//printf("Grad = (%lf, %lf).\n", gsl_vector_get(x, 0), gsl_vector_get(x, 1));

		if(status == GSL_SUCCESS && verbo)
			cout << "Achei um mínimo! Em: " << endl;
		if(verbo){
			cout << gsl_vector_get(ws->x, 0) << " " << gsl_vector_get(ws->x, 1) << " " << ws->f << endl;
			cout << endl;
		}

	}while(status == GSL_CONTINUE && iter < max_iter);

	if(status && verbo)
		cout << "Erro código " << status << ": " << gsl_strerror(status) << endl;

	if(iter == max_iter && verbo)
		cout << "Iterei muito e não achei... :(" << endl;

	if(verbo)
		cout << "No final devolvi " << ws->f << " no ponto " << gsl_vector_get(ws->x, 0) << " " << gsl_vector_get(ws->x, 1) << endl;

	double resultado;

	if(isnan(ws->f))
		resultado = 1e6;
	else
		resultado = ws->f;

	gsl_multimin_fdfminimizer_free(ws);

	return(resultado);
}

double simple_chisquare(int nBins, double N_pred[], double N_obs[]){
	double sigma2, result = 0;
	for(int i = 0; i < nBins; i++){
		if(N_obs[i] > 0)
			sigma2 = N_obs[i];
		else
			sigma2 = 1;
		result += pow(N_pred[i] - N_obs[i], 2)/sigma2;
	}
	return(result);
}