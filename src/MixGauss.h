#pragma once
#include "datastructure.h"

class MixGauss
{	public:
		MixGauss(void);
		~MixGauss(void);
	private:
		double *gausspara, *gaussprior, *lambda;
		int gausssz, lambdasz, groupn, *groupcode, *coderevmap;
		bool lessoneUpdate;
		const gsl_rng_type *T;
		gsl_rng *gammar;
		int maxxmsz, *xmsz, maxxvsz, maxymsz, *ymsz, maxyvsz, maxyxsz, **xmap, **ymap, *mapspace;

	public:
		int clustersz;
		double minerr, maxerr, tmpprop[10000];
		bool indc, nb;
		vector<vector<double> > modelparameter; 

	public:
		void computeLP(float *ydata, float *xdata, int id, double A, double *lp);

		void addPara(int g, float *yy, float *xx, int id);
		void removePara(int g, float *yy, float *xx, int id);
		void getStateCount(double *cn);
		void getStatePrior(double *&ppp, int &step, double A, double B);

		void updateLambda(double A);
		void initializePara(int clusterSZ, double *vh, float **dataYP, float **dataXP, int totalN, int L, int maxysz, int *ysz, int **ymp, int maxxsz, int *xsz, int **xmp, double A);
		void clearParameter();
		void rearrangeParameter(int *remap, int newclustersz, float **dataYP, float **dataXP, int totalN, int L, double A);
		void outputParameter(char *fout, double A);
		void updateParameterPrior(double A);
		
		void simData(int n, int id, int g, float *yp);

		double splitmergeCluster(int type, float **dataYP, float **dataXP, int totalN, int L, float *states, double A, int &mi, int &mj, int &tid);

		void imputeMissingTrack(int id, int L, float *yp, float *xp, float *yimp, float *ximp, float *states);//double *stateprob);

	private:
		double _dmvnorm(double *lambdap, float *yp, float *xp, int tymsz, int txmsz);
		void _getLambdaOne(double *ll, double *rr, double *tmpspace, double priorW, int i, int id);
		void _cholV(double const *A, int n, double *L);
		void _invL(double const *L, int n, double *iL, bool transpose);
		double _lgammaM(int q, double x);
		void _prepareGroupCode(int totalN);
		void _getMean(double const *pp, float *my, float *mx);
		void _gaussEM(double *pp, double *rr, double AA);
		void _collapsedPrediction(double const *pp, double const *phi, int id, double *rt, double *tmpspace, int &misyn, int *misy, int &misxn, int *misx, int *map);

		void _XtX(double const *X, int rn, int cn, double *X2);
		void _AtB(double const *A, double const *B, int ca, int cb, int n, double *C);
		void _ABt(double const *A, double const *B, int ra, int rb, int n, double *C);
		void _AB(double const *A, double const *B, int ra, int cb, int n, double *C);

		void _updateVh(double *Q, int K, double A, double *V);
};
