#ifndef SNPINFER_DATASTRUCTURE
#define SNPINFER_DATASTRUCTURE

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <vector>
#include <algorithm>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
using namespace std;

#define GSL_DLL
#include <gsl/gsl_sys.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_fit.h>


typedef struct SNPINFO {
        int chr, pos, label;
        int err, info, st, ed, shift;
	string snpid;
	double qual;
	vector<char> alleles;
} SNPINFO;

#define MINUSINFINITE -1000000000
#define PI 3.14159265359

typedef struct MYDATA {
	int indN, totalN, L, ploidity; //indN = N*plodity, totalN = N, but if replicates > 1, then totalN > N
	int *indIndex; //index of individuals (true individuals, not haploids, length+1 for total copies, =totalN) (necessary when there are replicates per individual)
	float *data; //length = indindex[last]*ploidity*L;
	int *asz;
	vector<vector<int> > fixState; //fix states, ~ pop, indN of them
	vector<vector<bool> > fixAllele; //fix data without imputation, ~ data, totalN of them; NOTE: false is fix, true is to impute!
	vector<vector<int> > breaks; //must contain start and end position as double breaks
	
	//additional info
	vector<SNPINFO> snpinfo;
	vector<string> indinfo;
	vector<string> fyinfo, fxinfo;
} MYDATA;

typedef struct TENSORPARA {
	//recombination parameters
	double priorW, probAlpha, defmut, *rPrior, *popr; //priorW is for recomb rate
							//probAlpha is for states proportion
	bool *rec;
	int *recombN, trueindN;

	//parameter for states
	int maxK, maxHapK; //maxHapK is used if number of states is fixed
	int *pop, *nh;	
	vector<vector<vector<float> > > param; //LxKxZ, where K=maxK, Z=asz
	vector<vector<int> > flips;
	double *Vh, *Q, A;

	//control parameters
	int burnin, mcmc;
	bool shrinkclass, singlecall, imputedata;
	bool addsample; //graudually add samples, used for haplotype inference
	bool samplefull; //true for haplotype inference
	char changeAlpha; //add (>0) or reduce (<0) probAlpha at the begining, add used for haplotype inference, reduce used for stratification mapping
	int splitstop; //# of pops, beyond which no more splitting
	int refweight; //weight of reference data with known states
	char maximum; //sample or maximize, 3bits abc: a: maximum recomb, b: maximum state, c: maximum allele
	bool samplemaximum; //indicating if first sample and then take maximum, used for pop inference
	bool logscale; //indicating if probability is in local scale, used for popinfer
	double recrate; //recombination rate to be used for set up recomb priors
	int recombcut; //above which position specific recomb is counted
	double heteroVh; //0 if Vh is global, > 0 if locally Vh changes. Only implemented for ploidity = 1, and not implemented in greedyswitch function
} TENSORPARA;


#endif
