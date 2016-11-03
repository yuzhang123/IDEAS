#include "genomicTensor.h"
#define STDFORMAT 1 

unsigned int timeA, timeB;
bool adjustB = true; //false: no adjust for B when update hyperparameter; true: adjust
bool useSoftState = false;//false; //use tmpSpace to calculate marginalized emission probability integrating out position state
bool stateinteraction = false;//individual state x position state, or independence

typedef struct MySortType {
        int index;
        double score, weight; 
} MySortType;

bool operator<(const MySortType &a, const MySortType &b)
{       return a.score < b.score;
} 

typedef struct MySortType2 {
        int id1, id2;
	string str;
} MySortType2;

bool operator<(const MySortType2 &a, const MySortType2 &b)
{       return a.str < b.str;
}
//things learned so far
//1) position specific Vh helps pops to converge
//2) taking imputation maximums helps gauss to converge
//3) hyper-parameters of priors should not change, although priors can be updated, otherwise priors may go to 0
//4) including position specific classes may need stronger priors, otherwise it will be hard to converge in small samples
//5) priors for emitMatrix should depend on distribution of gauss components, otherwise may be hard to converge in small samples

clock_t Time1, Time2;
double *TLP;

genomicTensor::genomicTensor(unsigned int rs):tensorHMMbase(rs)
{
timeA=timeB=0;
Time1=Time2=0;
TLP=NULL;
varstatus = NULL;

	dataYP = dataXP = NULL;
	dataRange = NULL;
	emitMatrix = NULL;
	gClass = NULL;
	posClassSZ = gClassSZ = 0;
	emitK = 0;
	maxG = 100;
	maxPosSZ = 100;
	independentInd = false;
	log2 = -1.;

	xmsz = ymsz = NULL;
	xmap = ymap = NULL;

	fixstateN = false;
}

genomicTensor::~genomicTensor()
{	if(dataYP != NULL) 
	{	for(int i = 0; i < mydata.totalN; i++) if(dataYP[i] != NULL) delete dataYP[i]; 
		delete dataYP;
	}
	if(dataXP != NULL) 
	{	for(int i = 0; i < mydata.totalN; i++) if(dataXP[i] != NULL) delete dataXP[i];
		delete dataXP;
	}
	if(dataRange != NULL) delete dataRange;
	if(emitMatrix != NULL) delete emitMatrix;
	if(gClass != NULL) delete gClass;

if(varstatus != NULL) delete varstatus;
	
	if(ymsz != NULL) delete ymsz;
	if(xmsz != NULL) delete xmsz;
	if(xmap != NULL) { for(int i = 0; i < mydata.totalN; i++) delete xmap[i]; delete xmap; }
	if(ymap != NULL) { for(int i = 0; i < mydata.totalN; i++) delete ymap[i]; delete ymap; }
}

void genomicTensor::_getAsz(MYDATA const &mydata)
{
	int i;
        for(i = 0; i < mydata.L; i++) mydata.asz[i] = mylik.clustersz;
}

/////////////////////////////////////////////////////////////////////////////////////
void genomicTensor::run(char const *finput, char const *fcov, int burnin, int mcmc, int maxHapK, int maxGG, int maxPos, double AA, double recr, double heteroVh, bool samplemaximum, char const *foutput, bool sqc, char const *fparam, bool add2, double merr, double mxerr, bool ind, bool indind, int startK, int fixC, vector<vector<int> > const &removelist, bool nb)
{

	int i, j;
	fixstateN = (bool)(fixC > 0);

	mydata.ploidity = 1;
	mydata.L = -1;
	initPara(tpara);
	tpara.priorW = 10.;
	tpara.maxHapK = maxHapK;
	if(maxGG > 0) maxG = maxGG;
	maxG = max(maxG, startK);
	if(maxPos > 0) maxPosSZ = maxPos;
	tpara.burnin = burnin;
	tpara.mcmc = mcmc;
	tpara.changeAlpha = 0;
	tpara.maximum = 0;//abc, recomb, state, allele
	tpara.samplemaximum = samplemaximum;//true;
	tpara.defmut = 0.1;
	tpara.A = 1.;
	tpara.recrate = recr;
	tpara.recombcut = 0;
	tpara.heteroVh = heteroVh; 
	tpara.samplefull = false; //fullsample is for haplotypes only, sample genome data by gClass

	mylik.minerr = merr;//min(0.01, overallsd / 10.);
	mylik.maxerr = mxerr;
	mylik.indc = ind;
	independentInd = indind;
	mylik.nb = nb;
	overallmean = overallmeanx = 0;
	overallsd = overallsdx = 1.;

	if(fparam != NULL) 
	{
		parseParameter(fparam);
		bool flag = false;
		for(i = (int)mylik.modelparameter.size() - 1; i > maxG; i--) 
			if((int)mylik.modelparameter[i].size() > 0)
			{	printf("cluster %d ignored as it is beyond the maximum # of clusters allowed (%d).\n", i-1, maxG);	
				flag = true;
			}
		if(flag) mylik.modelparameter.resize(maxG + 1);
	}

	readInput(mydata, finput, fcov, log2, add2); //need revision for fix state and fix alleles
	removeData(removelist);
//	printf("N=%d(%d), L=%d, p=%d\n", mydata.indN, mydata.totalN, mydata.L, ymsz);


//	if(mylik.maxerr == 100000000.) mylik.maxerr = max(1.,overallsd * 1.);
//	if(mylik.minerr > mylik.maxerr) mylik.minerr = mylik.maxerr;

/////////////////
	tpara.probAlpha = mydata.totalN/mydata.ploidity;//0.1;//min(1., max(1e-1, 10./ (double)mydata.totalN));
	if(AA > 0) tpara.probAlpha = AA;//tpara.probAlpha * AA; //changed 02/07/2015

//	tpara.A = tpara.probAlpha; //added on Jul 11th

	if(startK<=0) startK=20;
	tpara.A = (double)startK / 2.;
	if(fixstateN) { maxG = fixC; }
	mylik.clustersz = startK = min(startK, maxG); 
	posSZ = max(1,min(startK / 2, maxPosSZ));
	double tn[mylik.clustersz+1], tp[mylik.clustersz+1];
	for(i = 0; i < startK+1; i++) tn[i] = 0;
	_updateVh(tn, startK, max(tpara.A + (double)fixstateN * 1000000., (double)startK/2.), tp);
	mylik.initializePara(startK, tp, dataYP, dataXP, mydata.totalN, mydata.L, maxymsz, ymsz, ymap, maxxmsz, xmsz, xmap, tpara.priorW);

	double posRate = 0.1;
	basePop = new int[mydata.L];
	baseR = new bool[mydata.L];
	double tm[posSZ], tv[posSZ];
	for(i = 0; i < posSZ; i++) tm[i] = 0;
	_updateVh(tm, posSZ, tpara.A, tv);
	for(i = 1; i < posSZ; i++) tv[i] += tv[i - 1];
	for(i = 0; i < mydata.L; i++) 
	{	baseR[i] = false;
		if(i == 0 || gsl_ran_flat(gammar, 0, 1.) < posRate)
		{	double un = gsl_ran_flat(gammar, 0, tv[posSZ-1]);
			for(j = 0; j < posSZ; j++) if(un <= tv[j]) break;
			basePop[i] = j;
			baseR[i] = true;
		}
		else basePop[i] = basePop[i - 1];
	}

	tmpSpace = NULL;//new double[100 * mydata.L * maxG];
	inferStructure(mydata, tpara, foutput, sqc);
	if(tmpSpace != NULL) delete tmpSpace;
	
	outputResult(foutput, mydata);
printf("updatePrior1=%f, updateBase=%f\n", (double)timeA/1000000., (double)timeB/1000000.);
printf("baseForward=%f, baseBackward=%f\n", (double)Time1/1000000., (double)Time2/1000000.);
	
	if(mydata.data != NULL) delete mydata.data;
	if(mydata.indIndex != NULL) delete mydata.indIndex;
	
	delete basePop;
	delete baseR;
}

/////////////////////////////////////////////////////////////////////////////////////
void genomicTensor::_updatePara_remove(vector<float> &para, float data, float b, int id, int pos)
{	float *pdy = dataYP[id] + pos * ymsz[id], *pdx = NULL;
	if(maxxmsz > 0) pdx = dataXP[id] + pos * xmsz[id];
	para[(int)data] -= b;	
	mylik.removePara((int)data, pdy, pdx, id);
}

void genomicTensor::_updatePara_add(vector<float> &para, float data, float b, int id, int pos)
{	float *pdy = dataYP[id] + pos * ymsz[id], *pdx = NULL;  
	if(maxxmsz > 0) pdx = dataXP[id] + pos * xmsz[id];
	para[(int)data] += b;
	mylik.addPara((int)data, pdy, pdx, id);
}

//need to update pointers, bst and bed are breaks index, if they exist, otherwise they are positions
void genomicTensor::_sumforwardOne(int id, vector<int> const &breaks, int bst, int bed, TENSORPARA &tpara, int ploidity, float *&dp, float *&dp2, double const *&rp, int const *&nhp, vector<vector<vector<float> > >::const_iterator &ipm, int const *&aszp, double *&lpp, double *&lmm, double *&lss, int SSZ, int MK, vector<vector<int> > const &ss, vector<vector<int> > const &mlist, double *Vh2)
{	int i, j;
	double f, t, maxlp;	

double aaa = tpara.heteroVh;
double NNN = aaa * (double)mydata.indIndex[mydata.indN / mydata.ploidity] * mydata.ploidity + 1.;

	int ided = mydata.indIndex[id+1], idst = mydata.indIndex[id];
	TLP = NULL;
	if(*dp < 0)
	{	TLP = new double[(ided-idst) * *aszp];
		int tlpzz=0;
		float *tdp = dp;
		float *dpx = NULL, *dpy;
		for(i = idst; i < ided; i++, tdp += mydata.L, tlpzz += *aszp)
		{	dpy = dataYP[i] + bst * ymsz[i];
			if(maxxmsz > 0) dpx = dataXP[i] + bst * xmsz[i];
			mylik.computeLP(dpy, dpx, i, tpara.priorW, &TLP[tlpzz]);
		}
	}

	(*lss) = 0;
	if(bst <= 0)
	{	for(j = 0; j < SSZ; j++, lpp++)
        	{       (*lpp) = 0;
			f = _getLp(id, breaks, bst, false, ploidity, tpara, dp, dp2, &ss[j][0], &ss[j][1], ipm, nhp, aszp);
			//t = log(tpara.Vh[ss[j][0]]);
			int state = ss[j][0];
			t = log((tpara.Vh[state] + aaa*(double)(state < tpara.maxK)*(double)(*(nhp+state))) / NNN);
			(*lpp) = t + f;
			if(independentInd && j != id) { (*lpp) = MINUSINFINITE; }
			if(j == 0 || (*lpp) > maxlp) maxlp = (*lpp);
        	}
		lpp -= SSZ;
		for(j = 0; j < SSZ; j++, lpp++)
			(*lss) += exp((*lpp) - maxlp);
		(*lss) = log(*lss) + maxlp;
	}
	else
	{	double recp = *rp, nrecp = 1. - recp;
		double *tmp;
		if(tpara.singlecall) { recp = 1.-1e-10; nrecp = 1e-10; }//recp2 = 1.-1e-10;nrecp2 = 1e-10; }
		recp/=NNN;
		
	        double *tlpp = lpp - SSZ;
		if(TLP==NULL)
		{
			int i, j, asz = *aszp;
			double prop[2*(asz + 1)];

			double a = tpara.probAlpha * 1., ss = 1.;
			double *ppp;
			int step, sst = 0, o = *(basePop+bst);
			mylik.getStatePrior(ppp, step, tpara.A + (double)fixstateN * 1000000., tpara.priorW);
			if(emitMatrix != NULL) 
			{	if(useSoftState && tmpSpace != NULL)
				{	if(stateinteraction) 
					{	ppp = tmpSpace + ((emitK + 1) * bst) * (mylik.clustersz + 1);
						sst = mylik.clustersz + 1;
					}
					else ppp = tmpSpace + bst * (mylik.clustersz + 1);
					ss = ppp[mylik.clustersz];
				}
				else	
				{	if(stateinteraction) 
					{	ppp = emitMatrix + o * mylik.clustersz;
						sst = posSZ * mylik.clustersz;
					}
					else ppp = emitMatrix + (emitK * posSZ + o) * mylik.clustersz;
					ss = 1.;
				}
				step = 1;
			}
			ss*=a;
			double *statep = new double[SSZ * (mylik.clustersz + 1)];
			float *ddp = dp;
			
			for(i = idst; i < ided; i++, ddp += mydata.L)
			{	double *statepd = &statep[(int)(*ddp)], *statepa = &statep[mylik.clustersz], *pppd = ppp + (int)(*ddp) * step;
				for(j = 0; j < SSZ; j++, statepd += mylik.clustersz + 1, statepa += mylik.clustersz + 1)
				{	if(j < tpara.maxK && nhp[j]> 0)
					{	*statepd = (double)(*ipm)[j][(int)(*ddp)] + a*(*pppd);
						*statepa = nhp[j] + ss;
					}
					else
					{	*statepd = a*(*pppd);
						*statepa = ss;
					}
					if(j < SSZ - 1) pppd += sst;
				}
			}
			double sslp[SSZ];
			for(i = 0; i < SSZ; i++) sslp[i] = 1.;
			ddp = dp;
			for(i = idst; i < ided; i++, ddp += mydata.L)
			{	double *statepd = &statep[(int)(*ddp)], *statepa = &statep[mylik.clustersz], *pppd = ppp + (int)(*ddp) * step;
				for(j = 0; j < SSZ; j++, statepd += mylik.clustersz + 1, statepa += mylik.clustersz+1)
				{	sslp[j] *= *statepd / (*statepa);
					if(i < ided - 1)
					{	(*statepd) += 1;
						(*statepa) += 1;
					}
				}
			}
			delete statep;

			for(j = 0; j < SSZ; j++, lpp++, tlpp++)
			{	
        		        t = log((tpara.Vh[j] + aaa*(double)(j<tpara.maxK)*(double)(*(nhp+j))) * recp + exp((*tlpp) - (*(lss-1))) * nrecp) + (*(lss-1));
	        	        (*lpp) = t + log(sslp[j]);
				if(independentInd && j != id) { (*lpp) = MINUSINFINITE;  }
				if(j == 0 || maxlp < (*lpp)) maxlp = (*lpp);
			}
		}
		else
		{
			for(j = 0; j < SSZ; j++, lpp++, tlpp++)
			{	
				f = _getLp(id, breaks, bst, false, ploidity, tpara, dp, dp2, &ss[j][0], &ss[j][1], ipm, nhp, aszp);

  	      	        	int state = ss[j][0];
        		        t = log((tpara.Vh[state] + aaa*(double)(state<tpara.maxK)*(double)(*(nhp+state))) * recp + exp((*tlpp) - (*(lss-1))) * nrecp) + (*(lss-1));
//printf("<%d:%f>", j, f);
	        	        (*lpp) = t + f;
				if(independentInd && j != id) { (*lpp) = MINUSINFINITE;  }
				if(j == 0 || maxlp < (*lpp)) maxlp = (*lpp);
			}
//printf("\n");
		}
		tlpp = lpp - SSZ;
//		if(bst==74) printf("id=%d rec=%f (%d) ", id, recp, (int)tpara.pop[(1-id)*mydata.L+bst]); 
		for(j = 0; j < SSZ; j++, tlpp++)
		{	(*lss) += exp((*tlpp) - maxlp);
//			if(bst==74) printf("%f(%f),", *tlpp, tpara.Vh[j]);
		}
//		if(bst==74) 
//		{	for(j = idst; j < ided; j++) printf("%d,",(int)mydata.data[j*mydata.L+bst]);
//			printf("\n");
//		}
		(*lss) = log(*lss) + maxlp;
	}

	if(false && (bst==0 || bst == mydata.L - 1))
	{	double *tlpp = lpp-SSZ;
		printf("%d: ", id);
		for(i = 0; i < SSZ; i++, tlpp++) printf("%f(%f), ", *tlpp, log(tpara.Vh[i]));
		printf("| %f\n", *lss);
	}

	if(TLP != NULL) { delete TLP; TLP = NULL; }

	lmm += MK;
	lss++;
	ipm++;
	nhp += tpara.maxK;
	aszp++;
	rp++;
	dp++;
}

int genomicTensor::_tracebackOne(int id, vector<int> const &breaks, int bst, int bed, TENSORPARA &tpara, int ploidity, int *&pp, int *&pp2, bool *&rr, bool *&rr2, double const *&rp, int const *&nhp, vector<vector<vector<float> > >::const_iterator &ipm, double *&lpp, double *&lmm, double *&lss, int SSZ, int MK, vector<vector<int> > const &ss, vector<vector<int> > const &mlist, double *Vh2)
{
	pp--;
	nhp -= tpara.maxK;
	ipm--;
	lpp -= SSZ;
	lmm -= MK;
	lss --;
	rr--;
	rp--; 
	return 0;
}

double genomicTensor::_getLp(int id, vector<int> const &breaks, int bid, bool flip, int ploidity, TENSORPARA const &tpara, float const *dp, float const *dp2, int const *pp, int const *pp2, vector<vector<vector<float> > >::const_iterator ipm, int const *nhp, int const *aszp)
{
	int i, j, asz = *aszp;
	int ided = mydata.indIndex[id+1], idst = mydata.indIndex[id];
	double prop[2*(asz + 1)];

	_getProp(tpara, bid, (*pp), asz, ipm, nhp, basePop + bid, prop);
	////////////////////////
	
	double rt = 0, flp[asz], *tlp = TLP;
	for(i = idst; i < ided; i++, dp += mydata.L, tlp += asz)
	{	
if(TLP != NULL)
{		int k = 0;
		double maxlp;
		for(j = 0; j < asz; j++)
		{	
			flp[j] = tlp[j] + prop[j] - prop[asz];
			if(j == 0 || maxlp < flp[j]) { maxlp = flp[j]; k = j; }
		}
		double s = 0;
		for(j = 0; j < asz; j++)
			s += exp(flp[j] - maxlp);
		rt += log(s) + maxlp;
		if(i < ided - 1)
		{	prop[asz+1+k]+=1.;
			prop[k] = log(prop[asz+1+k]);
			prop[asz+1+asz]+=1.;
			prop[asz] = log(prop[asz+1+asz]);
		}	
	//	printf("%d: %f, ", (int)(*pp), rt);
	//	for(j = 0; j < asz; j++) printf("%f,",prop[j]-prop[asz]);
	//	printf("\n");
}
else
{ int k = (int)(*dp);
rt+=prop[k]-prop[asz];
if(i<ided-1)
{	prop[asz+1+k]+=1.;
	prop[k] = log(prop[asz+1+k]);
	prop[asz+1+asz]+=1.;
	prop[asz] = log(prop[asz+1+asz]);
}
}
	}

	return(rt);
}

double genomicTensor::_imputeOne(int id, vector<int> const &breaks, int bst, int bed, int flip, int ploidity, TENSORPARA const &tpara, float *dp, float *dp2, int const *pp, int const *pp2, vector<vector<vector<float> > >::const_iterator ipm, int const *nhp, int const *aszp)
{
	int i, j, asz = *aszp;
	int ided = mydata.indIndex[id+1], idst = mydata.indIndex[id];
/*
for(i=idst;i<ided;i++, tsp += asz)
{
		*(dp + (i - idst) * mydata.L) = (float)_sample(tsp, asz, tpara.maximum & 1);;
}
return;
*/

	double prop[2*(asz + 1)];
	_getProp(tpara, bst, (*pp), asz, ipm, nhp, basePop + bst, prop);
	////////////////////////////
	
	double a = tpara.priorW;
	double RT = 0;
	for(i = idst; i < ided; i++, dp += mydata.L)
	{	double tlp[asz], maxlp;
		float *dpy, *dpx = NULL;
		dpy = dataYP[i] + bst * ymsz[i];
		if(maxxmsz > 0) dpx = dataXP[i] + bst * xmsz[i];
int t = i * mydata.L + bst;
{
double un=0;//gsl_ran_flat(gammar, 0, 1.);
 {
		mylik.computeLP(dpy, dpx, i, a, tlp);
/*vector<double> TLP(asz),TLP1(asz);
for(j=0;j<asz;j++) TLP[j]=tlp[j],TLP1[j]=prop[j]-prop[asz]+tlp[j];
sort(TLP.begin(), TLP.end());
sort(TLP1.begin(),TLP1.end());
printf("<%f,%f> ", TLP[asz-1]-TLP[asz-2],TLP1[asz-1]-TLP1[asz-2]);
*/
//double ttlp[asz];for(j=0;j<asz;j++)ttlp[j]=tlp[j];
int oj=0;

		for(j = 0; j < asz; j++)
		{	tlp[j] += prop[j] - prop[asz];
			if(j == 0 || maxlp < tlp[j]) { maxlp = tlp[j]; oj = j; }
		}
	//	if((tpara.maximum & 1) == 0)
double sum = 0;
		{	for(j = 0; j < asz; j++)
			{	tlp[j] = exp(tlp[j] - maxlp);
sum += tlp[j];
			}
		}
int o = (int)(*dp);
double olp=tlp[o];tlp[o]=max(0.,tlp[o]-un);
		RT += maxlp;
		j = _sample(tlp, asz, tpara.maximum & 1);
if(false)//id==19 && bst == 176)//tpara.maximum > 0 || j >= asz)
{	printf("\n");
	for(t = 0; t<(asz+1);t++) printf("%d:%f, ",t,prop[t]);
	printf("| %d,%d\n", j,asz);
//	for(t = 0; t < asz; t++) printf("%d:%f ",t, ttlp[t]);
//	printf("\n");fflush(stdout);
/*	for(t = 0; t<(asz+1);t++) printf("%f,",prop[asz+1+t]);
	printf("\n");
	for(t = 0; t<asz;t++) printf("%f,",*(emitMatrix+(emitK*posSZ+basePop[bst])*mylik.clustersz+t));
	printf("\n");
	for(t = 0; t<asz;t++) printf("%f,",(*ipm)[*pp][t]);
	printf("\n");
	printf("%d:%d, state=%d\n", i, bst, (int)(*pp));
	exit(0);
*/
}
		*dp = (float)j;
		//*(dp + (i - idst) * mydata.L) = (float)j;
 }
}
		if(i < ided - 1)
		{	prop[asz+1+j]+=1.;
			prop[j] = log(/*exp(prop[j])*/prop[asz+1+j]);
			prop[asz+1+asz]+=1.;
			prop[asz] = log(/*exp(prop[asz])*/prop[asz+1+asz]);
		}	
	}
	return(RT);
}

void genomicTensor::_getProp(TENSORPARA const &tpara, int pos, int state, int asz, vector<vector<vector<float> > >::const_iterator ipm, int const *nhp, int const *op, double *prop)
{
	int i;
	
	////////////////////////////
	double a = tpara.probAlpha * 1., ss = 1.;
	double *ppp;
	int step;
	mylik.getStatePrior(ppp, step, tpara.A + (double)fixstateN * 1000000., tpara.priorW);

	if(emitMatrix != NULL) 
	{	if(useSoftState && tmpSpace != NULL)
		{	if(stateinteraction) ppp = tmpSpace + ((emitK + 1) * pos + min(state, emitK)) * (mylik.clustersz + 1);
			else ppp = tmpSpace + pos * (mylik.clustersz + 1);// + emitK/*min(tpara.maxK, state)*/ * mylik.clustersz;
			ss = ppp[mylik.clustersz];
		}
		else	
		{	if(stateinteraction) ppp = emitMatrix + (min(state, emitK) * posSZ + (*op)) * mylik.clustersz;
			else ppp = emitMatrix + (emitK * posSZ + (*op)) * mylik.clustersz;
			ss = 1.;
		}
		step = 1;
	}

	ss*=a;
	if(state < tpara.maxK && *(nhp + state)>0)
	{	vector<float>::const_iterator ipp = (*ipm)[state].begin();
		for(i = 0; i < asz; i++, ppp += step, ipp++) 
			prop[i] = ((*ipp) + a * (*ppp));
		prop[asz] = (double)*(nhp + state) + ss;
	}
	else
	{	
		for(i = 0; i < asz; i++, ppp += step) prop[i] = a * *ppp;
		prop[asz] = ss;
	}
//	printf("%d %d: ",pos, state);
	for(i = 0; i < asz + 1; i++) 
	{	prop[i + asz+1] = prop[i];
		prop[i] = log(prop[i]);
//		printf("%f,",prop[i]);
	}
//	printf("\n");
}

void genomicTensor::_updatePrior(vector<bool> const &init, bool oneround)
{	int i, j, k;


	if(oneround) 
	{	tpara.A = 1.;
	
		printf("cSZ:%d", mylik.clustersz);

//for(i = 0; i < 10; i++)
{
if(gID > tpara.burnin*1/4  && gID < tpara.burnin) 
{		_splitmergeCluster(0);
}
if(gID <= tpara.burnin *1/ 2)//gID >= tpara.burnin*3/5 && gID < tpara.burnin) 
{		_splitmergeCluster(1);
}
}
		if(gID <= tpara.burnin) _rearrangeState(1+1+(int)(tpara.A+0.5), maxG);
		printf("->%d ", mylik.clustersz);fflush(stdout);
	}

//printf("\n");
//for(i=0;i<mydata.L;i++) if(mydata.data[3*mydata.L+i]==4) printf("%d,",i);
//printf("\n");
//	double omaxerr = mylik.maxerr;
//	if(gID < tpara.burnin / 2) mylik.maxerr = mylik.minerr + (double)gID / (double)tpara.burnin * 2.;
	mylik.updateLambda(tpara.priorW);
//	mylik.maxerr = omaxerr;

	//update emission matrix	
	if(oneround)
	{	
/*
if(gID >= tpara.burnin)
{
float tmpy[mydata.L * (maxymsz - ymsz[0]) + 1], tmpx[mydata.L * (maxxmsz - xmsz[0]) + 1];
mylik.imputeMissingTrack(0, mydata.L, dataYP[0], dataXP[0], tmpy, tmpx, mydata.data);
float *py = &tmpy[0];
FILE *f = fopen("2.txt","a");
for(i=0;i<mydata.L;i++)
{	for(j=0;j<maxymsz-ymsz[0];j++,py++)
		fprintf(f, "%f ", *py);
	fprintf(f, "\n");
}
fclose(f);
}
*/
//if(gID < tpara.burnin) maxPosSZ = gID / 2 + 1; //this is to avoid huge memory usage at the begining iterations when data is too large
//else maxPosSZ = 100;
		if(gID < tpara.burnin && tpara.samplemaximum) tpara.maximum = 3 * (int)(gsl_ran_flat(gammar, 0, 1.) < min(0.5,(double)gID / tpara.burnin));	
//tpara.probAlpha = max(1, tpara.burnin - gID);
clock_t cst, ced;
cst=clock();

		int K = min(tpara.maxK, tpara.maxHapK);
		int pc2 = posSZ * mylik.clustersz, esz = (K + 1) * pc2;
		double *nemitMatrix = new double[esz];
		double esum[(K + 1)* posSZ];
		int n = mydata.indIndex[mydata.indN / mydata.ploidity];
		int idi = 0, idnn = n * mydata.ploidity;
		///////////////////////
		vector<vector<vector<float> > >::iterator ipm = tpara.param.begin();
		int *op = basePop, *nhp = tpara.nh;
		for(j = 0; j < esz; j++) nemitMatrix[j] = 0;
		for(j = 0; j < (K+1)*posSZ; j++) esum[j] = 0;	
		double *ep;

//double *logp = new double[(K+1)*posSZ*mylik.clustersz], *logpN = new double[(K+1)*posSZ];
//for(j = 0; j < (K+1)*posSZ*mylik.clustersz; j++) logp[j] = 0;
//for(j = 0; j < (K+1)*posSZ; j++) logpN[j] = 0;

		double *priorpp;
		int priorstep;
		mylik.getStatePrior(priorpp, priorstep, tpara.A + (double)fixstateN * 1000000., tpara.priorW); 

		for(j = 0; j < mydata.L; j++, ipm++, op++, nhp += tpara.maxK)
		{	ep = nemitMatrix + (*op) * mylik.clustersz;

			int tpop[idnn];
			int *tp = tpara.pop + j;
			int tttg[idnn];
			float *pdata = mydata.data + j;
			int ii = 0, jj = 1, l, xx = 0;
			for(k = 0; k < idnn; k++, pdata += mydata.L)
			{	if(k >= mydata.indIndex[jj]) { tp += mydata.L; jj++; }
				tpop[k] = *tp;
				tttg[k] = (int)(*pdata);
			}
			double PP[(K + 1) * mylik.clustersz];
			_getMLE_P(tpop, n, (*op), emitMatrix, emitK, K, PP, tttg, priorpp, priorstep);

			double sumk[mylik.clustersz], sumn = 0, nnn = 0;
			for(k = 0; k < mylik.clustersz; k++) sumk[k] = 0;
//double lpn=0;
			
			double *pPP = PP;
			for(k = 0; k < K; k++, ep += pc2)
				if(nhp[k] > 0)
				{	for(int l = 0; l < mylik.clustersz; l++, pPP++)
					{	ep[l] += (double)nhp[k]*(*pPP+1e-10);
						sumk[l] += (*ipm)[k][l];
//logp[k*pc2+(*op)*mylik.clustersz+l] += log(PP[k*mylik.clustersz+l] + 1e-100);
					}
//logpN[k*posSZ+(*op)]++;lpn++;
					esum[k*posSZ+(*op)]+=nhp[k];
					sumn += nhp[k];
					nnn++;
				}
				else pPP += mylik.clustersz;
			pPP=PP;
			for(int l = 0; l < K; l++)
				if(nhp[l] > 0)
				{	for(k = 0; k < mylik.clustersz; k++, pPP++)
					{	ep[k] += sumn * (*pPP+1e-10);
					}
				}
				else pPP+=mylik.clustersz;
		/*
			for(k = 0; k < mylik.clustersz; k++)
			{	for(int l = 0; l < K; l++) 
				{	ep[k] += sumn*(double)(nhp[l]>0)*(PP[l*mylik.clustersz+k]+1e-10);
//logp[K*pc2+(*op)*mylik.clustersz+k] += log(PP[l*mylik.clustersz+k] + 1e-100);
				}
			}	
		*/
//logpN[K*posSZ+(*op)]+=lpn;
			esum[K*posSZ+(*op)]+=sumn;//nnn;//1.;//sumn;
		}
/*
for(j = 0; j < posSZ; j++)
{	double ak[mylik.clustersz], M = tpara.probAlpha, sumak = 0, sumlp=0;
	for(k = 0; k < mylik.clustersz; k++) 
	{	ak[k] = nemitMatrix[K*pc2+j*mylik.clustersz+k]/(esum[K*posSZ+j]+1e-10);
//printf("%f(%f,%f);",ak[k], nemitMatrix[K*pc2+j*mylik.clustersz+k],esum[K*posSZ+j]);fflush(stdout);
		logp[K*pc2+j*mylik.clustersz+k] /= (logpN[K*posSZ+j]+1e-10);
		sumak += ak[k] * log(ak[k]);
		sumlp += ak[k] * logp[K*pc2+j*mylik.clustersz+k];
//printf("%f,",logp[K*pc2+j*mylik.clustersz+k]);
	}
//printf("\n");
	//_getMLE_HyperP1(&logp[K*pc2+j*mylik.clustersz], ak, M, mylik.clustersz, logpN[K*posSZ+j]);
	printf("<%d:%f> ", j, (double)(mylik.clustersz-1)*0.5772157/(sumak - sumlp));fflush(stdout);//this is maximum likelihood approximation of the precision parameter
//printf("\n");
//	printf("%d: ", j);
//	for(k = 0; k < mylik.clustersz; k++) nemitMatrix[K*pc2+j*mylik.clustersz+k] = ak[k]*esum[K*posSZ+j];//printf("%f(%f),", ak[k], nemitMatrix[K*pc2+j*mylik.clustersz+k]/esum[K*posSZ+j]);
//	printf("| %f\n", M);
}
delete logp; delete logpN;
*/

		ep = nemitMatrix;
		double totalep[mylik.clustersz], tttp[mylik.clustersz], tp[mylik.clustersz + 1];
		for(j = 0; j < mylik.clustersz; j++)
		{	totalep[j] = 0;
			for(k = 0; k < posSZ; k++) totalep[j] += ep[K*pc2+k*mylik.clustersz+j];
			totalep[j] *= tpara.probAlpha / ((double)adjustB*(double)(n - 1) + 1. + tpara.probAlpha);
		}
		_updateVh(totalep, mylik.clustersz, tpara.A, tttp, false);
		for(j = 0; j < (K+1)*posSZ; j++, ep += mylik.clustersz)
		{	for(k = 0; k < mylik.clustersz; k++) 
			{	ep[k] *= tpara.probAlpha / ((double)adjustB*(double)(n - 1) + 1. + tpara.probAlpha);
				ep[k] += 100.*tttp[k];
			}
			_updateVh(ep, mylik.clustersz, tpara.A, tp, false);
			for(k = 0; k < mylik.clustersz; k++) ep[k] = tp[k];
			esum[j] = 1.;
		}
		///////////////////////	

		if(emitMatrix != NULL) delete emitMatrix;
		emitMatrix = nemitMatrix;
		emitK = K;
ced=clock();
timeA+=(unsigned int)(ced-cst);
		//update basePop
cst=clock();
		_updateBase();
ced=clock();
timeB+=(unsigned int)(ced-cst);

		//sample state members
		if(false && gID >= tpara.burnin)
		{	int *pc;
			if(gClass == NULL) 
			{	gClassSZ = mylik.clustersz + 5;
				posClassSZ = posSZ + 5;
				unsigned int sz = gClassSZ * n * mydata.L + posClassSZ * mydata.L;
				gClass = new int[sz];
				pc = gClass;
				for(unsigned int j = 0; j < sz; j++, pc++) *pc = 0;
				maxG = min(maxG, gClassSZ);
				maxPosSZ = min(maxPosSZ, posClassSZ);
			}
			pc = gClass;
			float *pd = mydata.data;
			for(j = 0; j < n; j++)
			{	for(k = 0; k < mydata.L; k++, pc += gClassSZ, pd++)
					(*(pc + (int)(*pd))) = (*(pc + (int)(*pd))) + 1;
			}
			int *po = basePop;
			for(k = 0; k < mydata.L; k++, pc += posClassSZ, po++)
			{	(*(pc + *po)) = (*(pc + *po)) + 1;
			}
		}
	}

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void genomicTensor::readInput(MYDATA &mydata, char const *input, char const *fcov, double log2, bool add2)
{
	FILE *f;
	char tmp[100000];
	int i, j, k, n, l;
	
	_loadData(dataYP, input, mydata.indinfo, mydata.fyinfo, mydata.snpinfo, mydata.indN, mydata.L, mydata.totalN, maxymsz, ymsz, ymap, mydata.indIndex, overallmean, overallsd, log2);
	
	maxxmsz = 0; 
	if(fcov != NULL)
	{	vector<string> tmpindinfo;
		vector<SNPINFO> tmpsnpinfo;
		int N, L, tN, *tmpindIndex = NULL;
		_loadData(dataXP, fcov, tmpindinfo, mydata.fxinfo, tmpsnpinfo, N, L, tN, maxxmsz, xmsz, xmap, tmpindIndex, overallmeanx, overallsdx, log2);
		if(add2 && maxxmsz > 0)
		{	for(i = 0; i < mydata.totalN; i++)
			{	float *tdata = new float[L * xmsz[i] * 2];
				for(j = 0; j < xmsz[i]; j++)
				{	double m = 0;
					for(k = 0; k < L; k++) m += dataXP[i][k * xmsz[i] + j];
					m /= (double)L;
					for(k = 0; k < L; k++) 
					{	tdata[k * xmsz[i] * 2 + j * 2] = dataXP[i][k * xmsz[i] + j];
						tdata[k * xmsz[i] * 2 + j * 2 + 1] = pow((double)dataXP[i][k * xmsz[i] + j] - m, 2.);
					}
				}
				delete dataXP[i];
				dataXP[i] = tdata;
				int *tmap = new int[xmsz[i] * 2];
				for(j = 0; j < xmsz[i]; j++) 
				{	tmap[j * 2] = xmap[i][j] * 2;
					tmap[j * 2 + 1] = xmap[i][j] * 2 + 1;
				}
				delete xmap[i];
				xmap[i] = tmap;
				xmsz[i] *= 2;
			}
		}
		if(tmpindIndex != NULL) tmpindIndex = NULL;
	}
	else
	{	xmsz = new int[mydata.totalN];
		xmap = new int*[mydata.totalN];
		for(i = 0; i < mydata.totalN; i++)
		{	xmsz[i] = 0;
			xmap[i] = NULL;
		} 
		dataXP = new float*[mydata.totalN];
		for(i = 0; i < mydata.totalN; i++) dataXP[i] = NULL;
	}

	printf("ysz=%d, xsz = %d, overall ymean=%f, overall ysd=%f\n", maxymsz, maxxmsz, overallmean, overallsd);
	mydata.data = new float[mydata.totalN * mydata.L];
	float *dp = mydata.data;
	for(i = 0; i < mydata.totalN * mydata.L; i++) *dp++ = -1.;
	mydata.indN *= mydata.ploidity;//requires preset of ploidity
}

void genomicTensor::outputResult(char const *fname, MYDATA const &mydata)
{
	int i, j, k;
	char str[200];
	mylik.clustersz = 0;
	sprintf(str, "%s.g", fname);
	FILE *f = fopen(str, "w");
	for(i = 0; i < mydata.L; i++)
	{	fprintf(f, "%s chr%d %d ", mydata.snpinfo[i].snpid.c_str(), mydata.snpinfo[i].chr, mydata.snpinfo[i].pos);
		for(j = 0; j < mydata.totalN; j++) 
		{	//fprintf(f, "%d ", (int)mydata.data[j * mydata.L + i]);
			int mm = (int)mydata.data[j*mydata.L+i];
			if(gClass != NULL)
			{	int *gp = gClass + (j * mydata.L + i) * gClassSZ;
			 	for(k = 1; k < gClassSZ; k++) if(gp[mm] < gp[k]) mm = k;
				//fprintf(f, "%d(%d)%d ", mm, gp[mm], (int)mydata.data[j*mydata.L+i]);
			}
			fprintf(f, "%d ", mm);
			mydata.data[j * mydata.L + i] = mm;
			mylik.clustersz = max(mylik.clustersz, mm);
		}
		//fprintf(f, "%d\n", basePop[i]);
		int mm = (int)basePop[i];
		if(gClass != NULL)
		{	int *gp = gClass + mydata.totalN * gClassSZ * mydata.L + i * posClassSZ;
			for(k = 1; k < posClassSZ; k++) if(gp[mm] < gp[k]) mm = k;
			//fprintf(f, "%d(%d)%d\n", mm, gp[mm], basePop[i]);
		}
		fprintf(f, "%d\n", mm);
	}
	fclose(f);
	mylik.clustersz++;

	if(outputproportion && gClass != NULL)
	{	sprintf(str, "%s.gpp", fname);
		f = fopen(str, "w");
		for(i = 0; i < mydata.L; i++)
		{	fprintf(f, "%d ", i);
			for(j = 0; j < mydata.totalN; j++)
			{	int *gp = gClass + (j * mydata.L + i) * gClassSZ;
				for(k = 0; k < gClassSZ; k++) fprintf(f, "%d ", gp[k]);
			}
			fprintf(f, "\n");
		}
		fclose(f);
	}

	_refreshParameter(false);
	sprintf(str, "%s.para", fname);
	mylik.outputParameter(str, tpara.priorW);

for(int id = 0; id < mydata.totalN; id++)
{	if(maxymsz - ymsz[id] == 0) continue;
/*	double prop[2*(mylik.clustersz + 1)], maxlp;
	int const *pp = tpara.pop + id * mydata.L, *nhp = tpara.nh;
	vector<vector<vector<float> > >::const_iterator ipm = tpara.param.begin(); 
	float *dpy = dataYP[id], *dpx = dataXP[id];
	double *tlp = new double[mydata.L * mylik.clustersz], *ptlp = tlp;
	for(i = 0; i < mydata.L; i++, pp++, ipm++, nhp += tpara.maxK, dpy += ymsz[id], dpx += xmsz[id], ptlp += mylik.clustersz)
	{	
		_getProp(tpara, i, (*pp), mylik.clustersz, ipm, nhp, basePop + i, prop);
		mylik.computeLP(dpy, dpx, id, tpara.priorW, ptlp);
		for(j = 0; j < mylik.clustersz; j++)
		{	ptlp[j] += prop[j] - prop[mylik.clustersz];
			if(j == 0 || maxlp < ptlp[j]) { maxlp = ptlp[j]; }
		}
//		if(id==5) printf("%d: ", i);
		for(j = 0; j < mylik.clustersz; j++)
		{	ptlp[j] = exp(ptlp[j] - maxlp);
//			if(id==5) printf("%f,",ptlp[j]);
		}
//		if(id==5) printf("\n");
	}
*/
	float *tmpdata = new float[mydata.L * (maxymsz + maxxmsz - ymsz[id] - xmsz[id]) + 1], *tmpy = tmpdata, *tmpx = tmpy + mydata.L * (maxymsz - ymsz[id]);
	mylik.imputeMissingTrack(id, mydata.L, dataYP[id], dataXP[id], tmpy, tmpx, mydata.data + id * mydata.L);//tlp);
//	delete tlp;

	char str[100];
	sprintf(str, "%s_impute.%s", fname, mydata.indinfo[id].c_str());
	f = fopen(str,"w");
	for(i = 0; i < maxymsz; i++)
	{	for(j = 0; j < ymsz[id]; j++) if(i == ymap[id][j]) break;
		if(j >= ymsz[id]) fprintf(f, "%s ", mydata.fyinfo[i].c_str());
	}
	for(i = 0; i < maxxmsz; i++)
	{	for(j = 0; j < xmsz[id]; j++) if(i == xmap[id][j]) break;
		if(j >= xmsz[id]) fprintf(f, "%s ", mydata.fxinfo[i].c_str());
	}
	fprintf(f, "\n");
	float *py = tmpy, *px = tmpx;
	for(i = 0; i < mydata.L; i++)
	{	for(j = 0; j < maxymsz - ymsz[id]; j++, py++) fprintf(f, "%f ", *py);
		for(j = 0; j < maxxmsz - xmsz[id]; j++, px++) fprintf(f, "%f ", *px);
		fprintf(f, "\n");
	}
	fclose(f);
	delete tmpdata;
}

/*
	if(varstatus == NULL) varstatus = new double[mydata.L * 3];
	vector<int> response(mydata.indN);
	for(i = 0; i < mydata.indN; i++) response[i] = i;
	mylik.updateLambda(tpara.priorW * 1e-3);
	_testVariability(response, varstatus);
*/
	if(varstatus != NULL)
	{	sprintf(str, "%s.var", fname);
		f = fopen(str, "w");
		for(i = 0; i < mydata.L; i++)
		{	fprintf(f, "%f %f %f\n", varstatus[i], varstatus[i + mydata.L], varstatus[i + 2 * mydata.L]);	
		}
		fclose(f);
	}
}

double genomicTensor::_lgammaM(int q, double x)
{
	int i;
	double rt = (double)q * (double)(q - 1) / 4. * log(PI);
	/*for(i = 0; i < q; i++)
	{	rt += gsl_sf_lngamma(x);
		x -= 0.5;
	}*/
	double s = gsl_sf_lngamma(x);
	for(i = 0; i < q; i += 2)
	{	rt += s;
		if(i + 2 < q) s -= log(x - (double)i / 2. - 1.);
	}
	s = gsl_sf_lngamma(x - 0.5);
	for(i = 1; i < q; i += 2)
	{	rt += s;
		if(i + 2 < q) s -= log(x - (double)i / 2. - 1.);
	}
	return rt;
}

void genomicTensor::_updateBase()
{
	int pz1 = posSZ + 1;
        int i, j, k, l, p2 = pz1 * pz1;
        double mn[posSZ], mp[pz1], tn[p2], trans[p2], a = (tpara.probAlpha);//changed 05/26/2015 //1.;//tpara.probAlpha;//changed 02/07/2015
	double posRate[pz1], rr[posSZ], rn[posSZ];

        for(i = 0; i < posSZ; i++) mn[i] = rr[i] = rn[i] = 0;
        for(i = 0; i < p2; i++) tn[i] = 0;//((double)((int)(i/posSZ)==(i%posSZ)) + 1./(double)posSZ) / 2.;
        int *op = basePop;
	bool *rp = baseR;
	mn[*op]++; op++; rp++;
        for(i = 1; i < mydata.L; i++, op++, rp++)
       	{	mn[*op]++;
		if(*rp)
		{	tn[(*(op-1)) * pz1 + (*op)]++;
       		}
		if(*op == *(op-1) && !(*rp))
		{	rr[*(op-1)]++;
		}
		rn[*(op-1)]++;
	}
	_updateVh(mn, posSZ, tpara.A, mp, false);
	if(posSZ >= maxPosSZ) mp[posSZ] = 0;
	for(i = 0; i < posSZ; i++) 
	{	posRate[i] = 1. - (rr[i]+tpara.defmut*tpara.priorW * 1.) / ((double)rn[i] + tpara.priorW * 1.);
	}
	posRate[posSZ] = tpara.defmut;

	double *np = tn, *tp = trans;
	for(i = 0; i < pz1; i++, np += pz1, tp += pz1)
	{	_updateVh(np, posSZ, tpara.A, tp);
		if(posSZ >= maxPosSZ) tp[posSZ] = 0; 
	}

	int mytotaln = mydata.indIndex[mydata.indN / mydata.ploidity];	
	int nn = mytotaln + 1;
	int pc2 = posSZ * mylik.clustersz;
	double *A = new double[(emitK + 1) * pc2 * nn], *app = A;
	double *epp = emitMatrix;
	double *priorpp;
	int priorstep;
	mylik.getStatePrior(priorpp, priorstep, tpara.A, tpara.priorW); 
	for(i = 0; i < emitK + 1; i++)
	{	for(j = 0; j < pc2; j++, epp++)
		{	double f = a * priorpp[(j % mylik.clustersz) * priorstep];
			if(emitMatrix != NULL) f = a * (*epp);
			(*app) = 0;//gsl_sf_lngamma(f);
			app++;
			for(k = 1; k < nn; k++, app++)
			{	(*app) = *(app-1) + log((double)k - 1. + f);
			}
		}
	}
	double B[nn];
	B[0] = 0;//gsl_sf_lngamma(a);
	for(k = 1; k < nn; k++) B[k] = B[k-1] + log((double)k - 1. + a);

	double A0[mylik.clustersz * nn];
	for(i = 0; i < mylik.clustersz; i++) 
	{	double f = a * priorpp[i * priorstep];
		A0[i * nn] = 0;//gsl_sf_lngamma(f);
		for(j = 1; j < nn; j++)
			A0[i * nn + j] = A0[i * nn + j - 1] + log((double)j - 1. + f);
	}

	vector<int> region(1,0);
	j = 0;
	for(i = 1; i < mydata.L; i++) if(mydata.snpinfo[i].chr != mydata.snpinfo[j].chr) { region.push_back(i); j = i; }
	region.push_back(i);
	int maxl = 0;
	for(i = 1; i < (int)region.size(); i++) maxl = max(maxl, region[i] - region[i - 1]);

	double *E = new double[2*pz1 * maxl], *M = E + pz1 * maxl;
        double *P = new double[(pz1 + 1) * maxl];
        vector<vector<vector<float> > >::iterator ipm = tpara.param.begin();
        int *nhp = tpara.nh;
	int *aszp = mydata.asz;

	int tspstep = 1, maxo = 0;
	if(stateinteraction) tspstep = emitK + 1;
	if(useSoftState)
	{	if(tmpSpace != NULL) delete tmpSpace;
		tmpSpace = new double[tspstep * mydata.L * (mylik.clustersz + 1)];
	}

	for(int r = 1; r < (int)region.size(); r++)
	{	int tL = region[r] - region[r - 1];
		
	        //forward
clock_t cst, ced;
cst=clock();
		double *pE = E, *maxP = P + pz1 * maxl; 
	        for(i = 0; i < tL; i++, P += pz1, ipm++, nhp += tpara.maxK, maxP++, aszp++)
        	{       double prob[pz1];
	                if(i==0)
        	        {       for(j = 0; j < pz1; j++) prob[j] = (mp[j]);
                	}
	                else
	                {       double *tP = P - pz1;
				double mmm = *(maxP-1);
                	        for(j = 0; j < pz1; j++) prob[j] = 0;
	                        double *tr = trans;
        	                for(j = 0; j < pz1; j++, tP++)
                	        {       double pr = posRate[j], npr = 1. - posRate[j], ttp = exp(*tP-mmm), ttpr = ttp * pr;
					for(k = 0; k < pz1; k++, tr++)
                                	        prob[k] += ttpr * (*tr);
	                        	prob[j] += ttp * npr;
				}
                	}
			for(j = 0; j < pz1; j++) prob[j] = log(prob[j]);

	                int asz = *aszp, K = emitK, jz=0, jzstep = mylik.clustersz * nn, Kz = K*pc2*nn;
        	        double emit[pz1], maxe;
			int numbers[K * asz], kid[K * asz], astep[K * asz], kn = 0, kzstep = pc2 * nn, k0 = 0, kh[K];
			for(k = 0; k < K; k++)
				if(nhp[k] > 0)
				{	vector<float>::iterator ippp = (*ipm)[k].begin();
					for(l = 0; l < asz; l++, ippp++)
					{	if((*ippp) > 0.5) { numbers[kn] = (int)(*ippp); kid[kn] = k; astep[kn] = l * nn; kn++; }
					}
					kh[k0++] = nhp[k];
				}
			for(j = 0; j < pz1; j++, jz+=jzstep)
			{	double eee = 0, *aaa;
				aaa = A + Kz + jz;
				if(j == posSZ) aaa = A0;
				for(k = 0; k < kn; k++)
				{	if(stateinteraction && (k == 0 || kid[k] != kid[k-1])) aaa = A + kid[k] * kzstep + jz;
					eee += aaa[astep[k] + numbers[k]];
				}
				for(k = 0; k < k0; k++)	eee -= B[kh[k]];
				emit[j] = eee;
       	         		if(j == 0  || maxe < eee) maxe = eee;
			}

			for(j = 0; j < pz1; j++) 
			{	P[j] = prob[j] + emit[j];//exp(emit[j] - maxe);
				*pE++ = emit[j];
				if(j == 0 || *maxP < P[j]) *maxP = P[j];
			}
        	}

ced=clock();Time1+=ced-cst;cst=ced;

		double *tmpP, emitW[posSZ + 1], gW[mylik.clustersz];
		int tspstep = 1;
		if(useSoftState)
		{	tmpP = tmpSpace + tspstep * (region[r] - 1) * (mylik.clustersz + 1);

			for(i = 0; i < mylik.clustersz; i++)
			{	gW[i] = 0;
			}
			for(i = 0; i < posSZ; i++)
			{	emitW[i] = 0;
				for(j = 0; j < mylik.clustersz; j++) emitW[i] += emitMatrix[emitK * pc2 + i * mylik.clustersz + j] * exp(gW[j] / (double)maxymsz);
			}
			emitW[posSZ] = 1;
		}

	        //backward
		double *pM = M + pz1 * (tL - 1);
		int oo = -1;
		op = basePop + region[r] - 1;
		rp = baseR + region[r] - 1;
		for(i = tL - 1; i >= 0; i--, op--, rp--)
		{	P = P - pz1;
			maxP--;
			pE -= pz1;
			double ss = 0;
			double prob[pz1], *tr = trans, ttt[2];
			if(i == tL - 1)
			{	for(j = 0; j < pz1; j++) prob[j] = exp(P[j] - *maxP);
				if(useSoftState) { for(j = 0; j < pz1; j++) { pM[j] = P[j], P[j] = pE[j]; /*ss+= pM[j];*/ } }
				ss = *maxP;
			}
			else
			{	double ms;
				ms = ttt[0] = P[oo] + log(1. - posRate[oo]);
				int jj = oo;
				for(j = 0; j < pz1; j++, jj += pz1) 
				{	prob[j] = P[j] + log(posRate[j] * trans[jj]);
					if(ms < prob[j]) ms = prob[j];
				}
				for(j = 0; j < pz1; j++) prob[j] = exp(prob[j] - ms);
				ttt[0] = exp(ttt[0] - ms);
				prob[oo] += ttt[0];
				ttt[1] = prob[oo];
				if(useSoftState)
				{
					int jj=0;
					for(j = 0; j < pz1; j++,jj+=pz1)
					{	double s = 0, ms, f[pz1], pr = posRate[j], npr = 1. - posRate[j];
						for(k = 0; k < pz1; k++)
						{	f[k] = log(((double)(j==k)*npr+pr*trans[jj + k])) + P[pz1+k];
							if(k == 0 || ms < f[k]) ms = f[k];
						}
						for(k = 0; k < pz1; k++)
						{	s += exp(f[k] - ms);
						}
						if(s > 0 && s < 100000);
						else
						{	printf("i=%d: ", i);
							for(int k = 0; k < pz1; k++) printf("%f(%f),",P[pz1+k], log(trans[j*pz1+k]));printf(" max=%f\n", *(maxP+1));
							exit(0);
						}
						s = log(s);
						pM[j] = P[j] + s;

						P[j] = pE[j] + s;
						if(j == 0 || ss < pM[j]) ss = pM[j];
					}
				}
			}
			if(useSoftState)
			{
				double *tsp = tmpP;
				int kkkk=0;
				for(int kkk = emitK + 1 - tspstep; kkk < emitK + 1; kkk++,kkkk+=pc2)
				{	double sumtsp = 0, *otsp = tsp, sw = 0, swn = 0, f;
					for(j = 0; j < mylik.clustersz; j++, tsp++)
					{	*tsp = 0;
						int tkk=j;
						for(int k = 0; k < pz1; k++,tkk+=mylik.clustersz)
						{	f = exp(pM[k] - ss);
							if(k < posSZ) (*tsp) += f * emitMatrix[kkkk+tkk]; 
							else (*tsp) += f * priorpp[j * priorstep];
							sw += f * emitW[k];
							swn += f;
						}
						sumtsp += *tsp;
						if(*tsp > 0 && *tsp < 100000);
						else
						{	for(int k = 0; k < pz1; k++) printf("%f(%f), ", pM[k], log(emitMatrix[kkk*pc2+k*mylik.clustersz+j])); printf("| %f, tsp=%f\n", ss, *tsp);
							for(int k = 0; k < pz1; k++) printf("%d:%f, ", k, pM[k]);
							exit(0);
						}
					}
					sw = max(0.1, 0. + sw / swn);
					sumtsp /= sw;
					for(j = 0; j < mylik.clustersz; j++) otsp[j] /= sumtsp;
					*tsp++ = sw;
				}
				tmpP -= tspstep * (mylik.clustersz + 1);
			}

			int ooo = oo;
			oo = (*op) = _sample(prob, pz1, tpara.maximum & 2);
			(*rp) = true;
			if(oo == ooo)
			{	if(gsl_ran_flat(gammar, 0, ttt[1]) <= ttt[0]) (*rp) = false;
			}
			if(oo == posSZ)
			{	if(!(*rp) && i < tL - 1) (*op) = *(op+1);
				else if((tpara.maximum & 1) == 0)
                			(*op) = min(maxPosSZ - 1, posSZ + (int)(- log(gsl_ran_flat(gammar, 0, 1.)) / log((1. + tpara.A) / tpara.A)));
			}
			maxo = max(maxo, *op);
			pM -= pz1;
		}
Time2+=clock()-cst;
	}

	////collapse basePop, update posSZ;
	maxo += 1;
	vector<bool> on(maxo, false);
	op = basePop;
	for(i = 0; i < mydata.L; i++, op++) on[*op] = true;
	vector<int> omap(max(posSZ, maxo), -1);
	int newpossz = 0;
	for(i = 0; i < maxo; i++) 
		if(on[i]) omap[i] = newpossz++;
	op = basePop;
	for(i = 0; i < mydata.L; i++, op++) *op = omap[*op];
	double *nemitMatrix = new double[(emitK + 1) * newpossz * mylik.clustersz];
	double *ep = emitMatrix, *nep = nemitMatrix;
	for(i = 0; i < emitK + 1; i++)
	{	int l = 0;
		for(j = 0; j < posSZ; j++, ep += mylik.clustersz)
			if(on[j])
			{	for(k = 0; k < mylik.clustersz; k++, nep++) *nep = ep[k];
				l++;
			}
		for(j = l; j < newpossz; j++)
			for(k = 0; k < mylik.clustersz; k++, nep++) *nep = priorpp[k * priorstep];
	}
	delete emitMatrix;
	emitMatrix = nemitMatrix;
	int oz = posSZ;
	posSZ = newpossz;
	printf(" posSZ=%d->%d(%d) ", oz, posSZ, tpara.maximum & 1); 
	//for(i=0;i<pz1;i++) printf("%4.3f(%4.3f),", mp[i],posRate[i]); printf("\n");
	//if(posSZ == 1) exit(0);
	
	delete P;
	delete A;

	delete E;
}

void genomicTensor::_getMLE_P(int *tpop, int N, int o, double *E, int oK, int K, double *PP, int *tttg, double *priorpp, int priorstep)
{	
	int i, j;
	int pc2 = posSZ * mylik.clustersz;
	double a = tpara.priorW;

//double *tppp;
//mylik.getStatePrior(tppp, i, tpara.A, tpara.priorW);

	double *prop[K + 1], *ip = PP;
	int step[K + 1];
	double nn[K + 1];
	for(i = 0; i < K + 1; i++)
	{	step[i] = 1;
		if(oK == 0) { prop[i] = priorpp; step[i] = priorstep; }
		else prop[i] = E + oK * pc2 + o * mylik.clustersz;
		for(j = 0; j < mylik.clustersz; j++, ip++) 
		{	*ip = 1e-10/(double)mylik.clustersz+(double)(gID > tpara.burnin / 2) * sqrt(tpara.probAlpha) * (*(prop[i] + j * step[i]));//changed 05/26/2015 
								//(double)(gID*0 > tpara.burnin / 2) * tpara.probAlpha * (*(prop[i] + j * step[i]));//changed 02/07/2015
//if(gID > tpara.burnin/2)			*ip += 100. * (double)mydata.totalN * tppp[j];
		}
		nn[i] = 1e-10 + (double)(gID > tpara.burnin / 2) * sqrt(tpara.probAlpha);//changed 05/26/2015
				//(double)(gID*0 > tpara.burnin / 2) * tpara.probAlpha;//changed 02/07/2015
//if(gID > tpara.burnin/2)		nn[i] += 100. * (double)mydata.totalN;
	}
	
	for(i = 0; i < N; i++)
	{	int ddd = tttg[i];//not compatible with ploidity
		 *(PP+tpop[i]*mylik.clustersz+ddd) = *(PP+tpop[i]*mylik.clustersz+ddd) + 1;
                nn[tpop[i]]++;
		*(PP+K*mylik.clustersz+ddd) = *(PP+K*mylik.clustersz+ddd)+1;
		nn[K]++;
	}
	
	ip = PP;
	for(i = 0; i < K + 1; i++)
	{	for(j = 0; j < mylik.clustersz; j++, ip++) (*ip) /= nn[i];
	}
}

void genomicTensor::_getMLE_HyperP(double *logp, double *ak, int sz)
{
	int i, j;
	double oak[sz];
	for(i = 0; i < sz; i++) oak[i] = ak[i] * tpara.probAlpha;

	double psiC = gsl_sf_psi(tpara.probAlpha);
	double tsiC = gsl_sf_psi_1(tpara.probAlpha);
	bool a_inner;
	int cn = 0;
	do
	{	a_inner = true;
		double dF[sz], q[sz], b = 0, iq = 0;
		for(i = 0; i < sz; i++)
		{	dF[i] = (psiC - gsl_sf_psi(oak[i]) + logp[i]);
			q[i] = -gsl_sf_psi_1(oak[i]);
			b += dF[i] / q[i];
			iq += 1./q[i];
		}
		double s = 0;
		for(i = 0; i < sz; i++)
		{	ak[i] = max(0., oak[i] - (dF[i] - b / (1./tsiC+iq)) / q[i]);
			s += ak[i];
		}
		double e2 = 0;
		for(i = 0; i < sz; i++)
		{	ak[i] /= s / tpara.probAlpha;
			e2 += (oak[i]-ak[i])*(oak[i]-ak[i])/oak[i];
			oak[i] = ak[i];
		}
		cn++;
		if(e2 < 0.01 || cn >= 10) a_inner = false;
	} while(a_inner);
}

void genomicTensor::_rearrangeState(int add, int maxG)
{	int i, j, k;
	vector<char> sel(mylik.clustersz, 0);
//add=3;//added on Jul 11th

	int n = mydata.indIndex[mydata.indN / mydata.ploidity] * mydata.ploidity;
	float *dp = mydata.data;	
	for(i = 0; i < n; i++)
	{	for(j = 0; j < mydata.L; j++, dp++)
		{	
			sel[(int)(*dp)] = 1;
		}
	}
	for(i = 1; i < (int)mylik.modelparameter.size(); i++)
		if((int)mylik.modelparameter[i].size() > 0) sel[i-1] = 2;

	int map[mylik.clustersz], remap[mylik.clustersz], remain[mylik.clustersz], newclustersz = 0;
	j = k = 0;
	for(i = 0; i < mylik.clustersz; i++) 
	{	if(sel[i] == 2) { map[i] = i; remap[i] = i; newclustersz = max(newclustersz, i + 1); }
		else if(sel[i] == 1) { while(sel[j] == 2) j++;  map[i] = j; remap[j] = i; j++; newclustersz = max(newclustersz, j); }
		else { map[i] = -1; remain[k++] = i; }
	}
	k = 0;
	for(i = j; i < mylik.clustersz; i++) { if(sel[i] != 2) remap[i] = remain[k++]; }
	int trueclustersz = max(newclustersz, min(newclustersz + add, maxG));
	
	if(newclustersz < mylik.clustersz || trueclustersz != mylik.clustersz)
	{	vector<vector<float> > empp(tpara.maxK, vector<float>(trueclustersz, 0));
		vector<vector<vector<float> > > newparam(mydata.L, empp);
		//for(i = 0; i < mydata.L; i++) newparam[i] = empp;
		vector<vector<vector<float> > >::iterator inp = newparam.begin();
		vector<vector<vector<float> > >::iterator ip = tpara.param.begin();
		for(i = 0; i < mydata.L; i++, ip++, inp++)
		{	
		//	(*inp) = empp;//.resize(tpara.maxK, vector<float>(trueclustersz, 0));
			vector<vector<float> >::iterator inpp = (*inp).begin();
			vector<vector<float> >::iterator ipp = (*ip).begin();
			for(j = 0; j < tpara.maxK; j++, ipp++, inpp++)
			{	for(k = 0; k < newclustersz; k++)
				{	(*inpp)[k] = (*ipp)[remap[k]];
				}
				for(k; k < trueclustersz; k++)
					(*inpp)[k] = 0;
			}
			mydata.asz[i] = trueclustersz;
		}
		tpara.param.clear();
		tpara.param = newparam;

		dp = mydata.data;
		for(i = 0; i < n; i++)
			for(j = 0; j < mydata.L; j++, dp++)
			{	if((int)(*dp) >= mylik.clustersz) printf("!!!%d,%d,%d\n",(int)(*dp), i, j),exit(0);
				(*dp) = (float)map[(int)(*dp)];
			}
	
		int oclustersz = mylik.clustersz;
		mylik.rearrangeParameter(remap, trueclustersz, dataYP, dataXP, mydata.totalN, mydata.L, tpara.priorW);
		_getAsz(mydata);

		if(emitMatrix != NULL)
		{	n = (emitK + 1) * posSZ;
			double *ep = emitMatrix;
			double *newemitMatrix = new double[n * trueclustersz], *nep = newemitMatrix;
			double tppp = 1. / (1. + tpara.A);;
			for(i = 0; i < n; i++, ep += oclustersz, nep += trueclustersz)
			{	double msum = 0;
				for(j = 0; j < trueclustersz; j++) 
				{	if(j < oclustersz) { nep[j] = ep[remap[j]]; msum += nep[j]; }
					else nep[j] = max(1e-100,(1. - msum)) * (1. - tppp) * pow(tppp, (double)(j - oclustersz));
				}
			}
			delete emitMatrix;
			emitMatrix = newemitMatrix;
		}
	}
}

void genomicTensor::_loadData(float **&datamatrix, char const *input, vector<string> &sid, vector<string> &fid, vector<SNPINFO> &snpinfo, int &indN, int &L, int &totalN, int &maxmsz, int *&msz, int **&map, int *&indIndex, double &mean, double &sd, double log2)
{
	ifstream f(input);
	string tmp;
	int i, j, k, n, l;
	
	if(f.bad()) printf("!!!cannot open file %s\n", input), exit(0);
	for(L = 0; getline(f, tmp); L++) ;
	L--;
	f.close();

	f.open(input);
	getline(f, tmp);
	l = (int)tmp.size(); l--; while(l>0 && (int)tmp[l] < 32) l--; l++;
	while(l > 0 && (tmp[l-1] == ' ' || tmp[l-1] == '\t')) l--;
	for(i = 0; i < l; i++) if(tmp[i] == ' ' || tmp[i] == '\t') break;
	for(i = i + 1; i < l; i++) if(tmp[i] == ' ' || tmp[i] == '\t') break;
	vector<string> sampleid, factorid, ids;
	while(i < l)
	{	for(j = i + 1; j < l; j++) if(tmp[j] == ' ' || tmp[j] == '\t') break;
		tmp[j] = 0;
		ids.push_back(&tmp[i + 1]);
		for(k = i + 1; k < j; k++) if(tmp[k] == '.') break;
		if(k < j) tmp[k] = 0;
		sampleid.push_back(&tmp[i + 1]);
		factorid.push_back(&tmp[k + (int)(k < j)]);
		if(k < j) tmp[k] = '.';
		if(j < l) tmp[j] = ' ';
		i = j;
	}
	n = (int)sampleid.size();
	sid.clear(); fid.clear();
	for(i = 0; i < n; i++)
	{	for(j = 0; j < (int)sid.size(); j++) if(sid[j] == sampleid[i]) break;
		if(j >= (int)sid.size()) sid.push_back(sampleid[i]);
	}	
	for(i = 0; i < (int)factorid.size(); i++)
	{	for(j = 0; j < (int)fid.size(); j++) if(fid[j] == factorid[i]) break;
		if(j >= (int)fid.size()) fid.push_back(factorid[i]);
	}	
	indN = (int)sid.size();
	maxmsz = (int)fid.size();
	vector<vector<int> > repn(indN, vector<int>(maxmsz, 0)), trepn = repn;
	for(i = 0; i < n; i++)
	{	for(j = 0; j < (int)sid.size(); j++) if(sid[j] == sampleid[i]) break;
		for(k = 0; k < (int)fid.size(); k++) if(fid[k] == factorid[i]) break;
		repn[j][k]++;
	}
	indIndex = new int[indN + 1];
	indIndex[0] = 0;
	for(i = 0; i < indN; i++)
	{	k = 0;
		for(j = 0; j < maxmsz; j++) k = max(k, repn[i][j]);
		indIndex[i + 1] = indIndex[i] + k;
	}
	
	totalN = indIndex[indN];
	msz = new int[totalN];
	map = new int*[totalN];
	for(i = 0; i < indN; i++)
	{	for(j = indIndex[i]; j < indIndex[i + 1]; j++)
		{	map[j] = new int[maxmsz];
			k = msz[j] = 0;
			for(l = 0; l < maxmsz; l++) 
				if(repn[i][l] > j - indIndex[i]) 
				{	msz[j]++;
					map[j][k++] = l;
				}
		}
	}
	
//	vector<int> tn(totalN, 0);
	vector<vector<int> > samplemap;
	for(i = 0; i < n; i++)
	{	int ss, ff, tt;
		for(ss = 0; ss < indN; ss++) if(sampleid[i] == sid[ss]) break;
		for(ff = 0; ff < maxmsz; ff++) if(factorid[i] == fid[ff]) break;
		vector<int> row(2, 0);
		row[0] = indIndex[ss] + trepn[ss][ff];
		for(tt = 0; tt < msz[row[0]]; tt++) if(map[row[0]][tt] == ff) break;
		row[1] = tt;//tn[indIndex[ss] + trepn[ss][ff]];
		samplemap.push_back(row); 
//		tn[indIndex[ss] + trepn[ss][ff]]++;
		trepn[ss][ff]++;
	}
	vector<string> newsid(totalN);
	for(i = 0; i < indN; i++)
	{	if(indIndex[i + 1] - indIndex[i] > 1)
		{	char tstr[1000];
			for(j = indIndex[i]; j < indIndex[i + 1]; j++)
			{	sprintf(tstr, "%s.%d", sid[i].c_str(), j - indIndex[i]); 
				newsid[j] = tstr;
			}
		}
		else newsid[indIndex[i]] = sid[i];
	}
	sid = newsid;
/*
for(i = 0; i < (int)samplemap.size(); i++)
{	printf("%d: %d,%d\n", i, samplemap[i][0], samplemap[i][1]);
}
printf("maxymsz=%d\n", maxmsz);
for(i=0;i<(int)fid.size();i++) printf("%s ", fid[i].c_str());printf("\n");
for(i=0;i<totalN;i++)
{	printf("%d: %d ", i, msz[i]);
	for(j = 0; j < msz[i]; j++) printf("%d,",map[i][j]);
	printf("\n");
}
for(i=0;i<indN;i++)
{	printf("%d: ", i);
	for(j=0;j<maxmsz;j++) printf("%d,",repn[i][j]);
	printf("\n");
}
*/
	snpinfo.clear();
	mean = sd = 0;
	int datan = 0;
	float row[n];
	datamatrix = new float*[totalN];
	float *dp[totalN];
	for(i = 0; i < totalN; i++) dp[i] = datamatrix[i] = new float[L * msz[i]];
	for(i = 0; getline(f, tmp); i++)
	{	l = (int)tmp.size();
		SNPINFO tsnp;
		char nd[100];
		sprintf(nd, "%d", i);
		tsnp.snpid = nd;
		if(tmp[0]=='c')
		{	tsnp.chr = atoi(&tmp[3]);
			if(tmp[3]=='X') tsnp.chr=23;
			else if(tmp[3]=='Y') tsnp.chr=24;
			else if(tmp[3] == 'M') tsnp.chr=25;
		}
		else tsnp.chr = atoi(&tmp[0]);
		for(j = 1; j < l; j++) if(tmp[j] == ' ' || tmp[j] == '\t') break;
		tsnp.pos = atoi(&tmp[j+1]);
		tsnp.qual = -1;
		tsnp.label = 1;
		snpinfo.push_back(tsnp);
		for(j = j + 1; j < l; j++) if(tmp[j] == ' ' || tmp[j] == '\t') break;
				
		k = 0;
		while(j < l - 1)
		{	double value = atof(&tmp[j + 1]);
			if(log2 > 0) value = log(value + log2) / log(2.);
			row[k++] = value; 
			for(j = j + 1; j < l - 1; j++) if(tmp[j] == ' ' || tmp[j] == '\t') break;
		}
		for(j = 0; j < n; j++) 
		{	mean += (double)row[j]; sd += (double)(row[j] * row[j]); 
			dp[samplemap[j][0]][samplemap[j][1]] = row[j];
		}
		datan += n;
		for(j = 0; j < totalN; j++) dp[j] += msz[j];
	}
	f.close();
	mean /= ((double)datan + 1e-10);
	sd = sqrt(max(1e-5, sd / ((double)datan + 1e-10) - mean * mean));
}

//solve psi(x)=y
double genomicTensor::_solveDigamma(double x0, double y, double precision, int maxitern)
{	int i = 0;
	double x = x0, err;
	while(i < maxitern)
	{	x = x0 - (gsl_sf_psi(x0) - y) / gsl_sf_psi_1(x0);
		err = fabs((x - x0) / x0);
		if(err < precision) break;
		x0 = x;
		i++;
	}
	return(x);
}

//mle of mean and precision parameters
void genomicTensor::_getMLE_HyperP1(double *logp, double *ak, double &M, int sz, double N, double precision, int maxitern)
{	int i, cn = 0;
	double newak[sz], newM;
	if(N<=0) return;
	while(cn < maxitern)
	{	double s, err, err2, y = 0;
		for(i = 0; i < sz; i++)
		{	
			y += ak[i] * (logp[i] - gsl_sf_psi(M * ak[i]));
		}
		s = 0;
		for(i = 0; i < sz; i++)
		{	
			newak[i] = _solveDigamma(ak[i], logp[i] - y, precision, maxitern);
			s += newak[i];
		}
		err = y = 0;
		for(i = 0; i < sz; i++) 
		{	newak[i] /= s;
			err += (ak[i] - newak[i]) * (ak[i] - newak[i]) / ak[i];
			ak[i] = newak[i];
			y += ak[i] * (logp[i] - gsl_sf_psi(M * ak[i]));
		}
		double ds = N * (gsl_sf_psi(M) + y);
		y = 0;
		for(i = 0; i < sz; i++)
			y += ak[i] * ak[i] * gsl_sf_psi_1(M * ak[i]);	
		double ds2 = N * (gsl_sf_psi_1(M) - y);
		newM = M;//min(100.,max(1.,M - ds / ds2));
		err2 = (M - newM) * (M - newM) / M;
		M = newM;
		if(err < precision && err2 < precision) break;
		cn++;
	}
}

void genomicTensor::parseParameter(char const *parafile)
{
	FILE *f = fopen(parafile, "r");
	if(f != NULL)
	{	char tmp[1000];
		while(fgets(tmp, 1000, f) != NULL)
		{	int i, j, k, l = (int)strlen(tmp) - 1;
			double f;
			if(l < 3) continue;
			string str(tmp);
			for(i = 0; i < l; i++) if(tmp[i] == '=') break;
			if(str.substr(0,6) == "burnin")
			{	j = atoi(&tmp[i + 1]);
				if(j > 0) tpara.burnin = j; 
				printf("burnin = %d\n", tpara.burnin);
			}
			else if(str.substr(0,4) == "mcmc")
			{	j = atoi(&tmp[i+1]);
				if(j > 0) tpara.mcmc = j;
				printf("mcmc = %d\n", tpara.mcmc);
			}
			else if(str.substr(0,13) == "samplemaximum")
			{	tpara.samplemaximum = (bool)atoi(&tmp[i+1]);
				printf("samplemaximum = %d\n", (int)tpara.samplemaximum);
			}
			else if(str.substr(0,6) == "maxN_S")
			{	tpara.maxHapK = atoi(&tmp[i+1]);
				printf("maxHapK = %d\n", tpara.maxHapK);
			}
			else if(str.substr(0,6) == "maxN_O")
			{	j = atoi(&tmp[i + 1]);
				if(j > 0) maxPosSZ = j;
				printf("maxPosSZ = %d\n", maxPosSZ);
			}
			else if(str.substr(0,6) == "maxN_K")
			{	j = atoi(&tmp[i+1]);
				if(j > 0) maxG = j;
				printf("maxG = %d\n", maxG);
			}
			else if(str.substr(0,4) == "log2")
			{	log2 = atof(&tmp[i+1]);
				printf("log2 = %f\n", log2);
			}
			else if(str.substr(0,1) == "A")
			{	f = atof(&tmp[i+1]);
				if(f > 0) tpara.A = f;
				printf("A = %f\n", tpara.A);
			}
			else if(str.substr(0,7) == "recrate")
			{	f = atof(&tmp[i+1]);
				if(f >= 0) tpara.recrate = f;
				printf("recrate = %f\n", tpara.recrate);
			}
			else if(str.substr(0,5) == "minsd")
			{	f = atof(&tmp[i+1]);
				if(f >= 0) mylik.minerr = f;
				printf("minerr = %f\n", mylik.minerr);
			}
			else if(str.substr(0,7) == "fixcoef")
			{	j = atoi(&tmp[i + 1]);
				for(k = i + 2; k < l; k++)
				{	if(tmp[k] == ' ' || tmp[k] == '\t') break;
				}
				vector<double> row;
				while(k < l)
				{	f = atof(&tmp[k + 1]);
					for(k = k + 1; k < l; k++) if(tmp[k] == ':') break;
					double w = atof(&tmp[k + 1]);
					for(k = k + 1; k < l; k++) if(tmp[k] == ' ' || tmp[k] == '\t') break;
					row.push_back(f);
					row.push_back(w);
				}
				if((int)mylik.modelparameter.size() <= max(-1,j)+1) mylik.modelparameter.resize(max(-1,j)+2, vector<double>());
				mylik.modelparameter[max(-1,j)+1] = row;
				printf("%d(%d): ", j, (int)mylik.modelparameter.size());
				for(k = 0; k < (int)row.size(); k+=2) printf("%f:%f ", row[k], row[k+1]);
				printf("\n");
			}
			
		}
		fclose(f);
	}
}

void genomicTensor::_splitmergeCluster(int type)
{	
	double *oemt = emitMatrix;
	emitMatrix = NULL;
//	double lp0 = _logP(mydata, tpara);
//	float *newdata = new float[mydata.totalN * mydata.L];
	int i, j, k, mi, mj, tid;
//for(i=0;i<mydata.totalN;i++) printf("%d ",(int)mydata.data[i*mydata.L+166]);printf("\n");
	double dlp = mylik.splitmergeCluster(type, dataYP, dataXP, mydata.totalN, mydata.L, mydata.data, tpara.priorW, mi, mj, tid);
//	double odlp = dlp;
	
	if(mi>=0)
	{
		if(type == 1) mj = tid;
		int estep = 0;
		double tepp[mylik.clustersz], *epp = &tepp[0];
		if(emitMatrix != NULL) 
		{	epp = emitMatrix + emitK * posSZ * mylik.clustersz;
			estep = mylik.clustersz;
		}
		else
		{	mylik.getStateCount(&tepp[0]);
			double ss = 1e-100;
			for(i = 0; i < mylik.clustersz; i++) ss += tepp[i];
			for(i = 0; i < mylik.clustersz; i++) tepp[i] /= ss;
		}
		int *op = basePop;
		double a = tpara.probAlpha * 1.;
		vector<vector<vector<float > > >::iterator ipm = tpara.param.begin();
		double N1 = 1., N2 = 1.;//, sign = 2. * (double)(type == 0) - 1.;
		for(i = 0; i < mydata.L; i++, ipm++, op++)
		{	double *eee = epp + estep * (*op);
			for(j = 0; j < tpara.maxK; j++)
			{	int n1 = (int)(*ipm)[j][mi], n2 = (int)(*ipm)[j][mj];
				if(n1 > 0 && n2 > 0)
				{	
					dlp += gsl_sf_lngamma(n1 + n2 + a * (*(eee + mi) + *(eee + mj)));
					dlp -= (gsl_sf_lngamma(n1 + a * *(eee + mi)) + gsl_sf_lngamma(n2 + a * *(eee + mj)));
					dlp -= gsl_sf_lngamma(a * (*(eee + mi) + *(eee + mj)));
					dlp += (gsl_sf_lngamma(a * *(eee + mi)) + gsl_sf_lngamma(a * *(eee + mj)));
				}
				N1+=n1;N2+=n2;
			}
		}
//dlp += gsl_sf_lngamma(N1+N2)-gsl_sf_lngamma(N1)-gsl_sf_lngamma(N2);

		double un = gsl_ran_flat(gammar, 0, 1.);

//printf("[ %f < %f ? %d]\n", log(un/(1.-un)), dlp, type);
		if(type == 0)
		{	
			if(log(un / (1. - un)) <= dlp) //to merge
			{	/*ipm = tpara.param.begin();
				for(i = 0; i < mydata.L; i++, ipm++)
				{	for(j = 0; j < tpara.maxK; j++)
					{	(*ipm)[j][tid] = (*ipm)[j][mi] + (*ipm)[j][mj];
						(*ipm)[j][mi] = (*ipm)[j][mj] = 0;
					}
				}*/
				float *dp = mydata.data;
				for(i = 0; i < mydata.totalN; i++)
					for(j = 0; j < mydata.L; j++, dp++)
						if((*dp) == mi || (*dp) == mj) (*dp) = tid;
			}
			_refreshParameter();
		}
		else
		{	if(log(un / (1. - un)) <= dlp) //fail to split
			{	float *dp = mydata.data;
				for(i = 0; i < mydata.totalN; i++)
					for(j = 0; j < mydata.L; j++, dp++)
						if((*dp) == tid) (*dp) = mi;
			}
			_refreshParameter();
		}
	}
	emitMatrix = oemt;
}

void genomicTensor::_getLpVar(float *ypp, float *xpp, int nst, int ned, int *tymsz, int *txmsz, vector<vector<double> > const &props, double priorW, int asz, vector<double> &rt)
{
	int i, j, k, psz = (int)props.size();
	int n = ned - nst;

	TLP = new double[n * asz];
	int tlpzz=0;
	double *dpx = NULL;
	for(i = 0; i < n; i++, ypp += tymsz[i], xpp += txmsz[i], tlpzz += asz)
		mylik.computeLP(ypp, xpp, nst + i, priorW, &TLP[tlpzz]);

	rt.clear();
	rt.resize(psz, 0);
	double *tlp = TLP, maxlp, flp[asz];
	for(i = 0; i < n; i++, tlp += asz)
	{	for(k = 0; k < psz; k++)
		{	for(j = 0; j < asz; j++)
			{	flp[j] = tlp[j] + props[k][j];
				if(j == 0 || maxlp < flp[j]) { maxlp = flp[j]; }
			}
			double s = 0;
			for(j = 0; j < asz; j++)
				s += exp(flp[j] - maxlp);
			rt[k] += log(s) + maxlp;
		}
	}

	delete TLP; TLP = NULL;
}

void genomicTensor::_testVariability(vector<int> const &response, double *vprob)
{       int i, j, k;
	vector<MySortType> groupid(mydata.indN);
	for(i = 0; i < mydata.indN; i++)
	{	groupid[i].index = i;
		groupid[i].score = (double)response[i];
	}
	sort(groupid.begin(), groupid.end());
	int rN = 0;
	for(i = 0; i < mydata.indN; i++)
		if(i == mydata.indN - 1 || groupid[i].score != groupid[i+1].score) rN++;
	
	vector<MySortType2> ggg;
	for(j = 0; j < mydata.L; j++)
	{	
		vector<vector<int> > tggg(rN);
		float *dp = mydata.data + j;
		int rid = 0;
		for(int ii = 0; ii < mydata.indN; ii++)
		{	i = groupid[ii].index;
			int idst = mydata.indIndex[i], ided = mydata.indIndex[i + 1];
			for(k = 0; k < ided - idst; k++, dp += mydata.L)
			{	tggg[rid].push_back((int)(*dp));	
			}
			if(ii == mydata.indN - 1 || groupid[ii].score != groupid[ii + 1].score) rid++;
		}
		for(i = 0; i < rN; i++)
		{	
			sort(tggg[i].begin(), tggg[i].end());
			MySortType2 my;
			my.id1 = j; my.id2 = i;
			for(k = 0; k < (int)tggg[i].size(); k++)
			{	//ostringstream convert;
				char tmpstr[100];
				//convert<<tggg[i][k];
				sprintf(tmpstr, "%d,", tggg[i][k]);
				my.str += tmpstr;//convert.str() + ",";
			}
			ggg.push_back(my);
		}
	}

	vector<vector<int> > glabel(mydata.L, vector<int>(rN, -1));
	sort(ggg.begin(), ggg.end());
	vector<double> nnn;
	j = 0;
	int ok = 0;
	int pN = (int)ggg.size();
	double psuedo = 0, cn = psuedo+1, totalcn = 0;//, cutvalue = 0;
	for(k = 1; k < (int)ggg.size(); k++) 
	{	if(ggg[k].str != ggg[j].str || k == (int)ggg.size() - 1)
		{	if(true)//cn > cutvalue)
			{	nnn.push_back(cn + (int)(k==pN-1));
				totalcn += cn;
				for(int l = ok; l < k + (int)(k==pN-1); l++) glabel[ggg[l].id1][ggg[l].id2] = j;
				j++;
			}
	//		else printf("cut=%d\n", (int)cn);
			ggg[j] = ggg[k];
			cn = psuedo;
			ok = k;
		}
		cn++;
	}
	ggg.resize(j);
	pN = j;
	
	vector<vector<double> > prop0(pN, vector<double>(mylik.clustersz, 0));
	for(i = 0; i < pN; i++)
	{	char str[100];
		sprintf(str, "%s", ggg[i].str.c_str());
		k = (int)strlen(str) - 1;
		j = 0;
		int c = 0;
		do
		{	prop0[i][atoi(&str[j])]++;	
			c++;
			for(j; j < k; j++) if(str[j] == ',') break;
			j++;
		} while(j < k);
//		for(j = 0; j < mylik.clustersz; j++) prop0[i][j] = log((prop0[i][j] + 1e-100) / (double)c);
	}

vector<vector<double> > score;
_pairEnrich(prop0, nnn, score);
printf("\n{\n");
for(i=0;i<(int)score.size();i++)
{	for(j=0;j<(int)score[i].size();j++)
		printf("%3.3f ", score[i][j]);
	printf("\n");
}
printf("}\n");
for(i = 0; i < pN; i++)
{	printf("%d: ", (int)nnn[i]);
	for(j = 0; j < mylik.clustersz; j++) printf("%d ", (int)prop0[i][j]);
	printf("\n");
}

	printf("pN=%d(%d): ", (int)nnn.size(), (int)prop0.size());for(i=0;i<(int)nnn.size();i++) printf("%d,",(int)nnn[i]);printf("\n");
	for(i = 0; i < pN; i++)
	{	printf("%d: ", i);
		for(j = 0; j < mylik.clustersz; j++) printf("%d,",(int)prop0[i][j]); printf("\n");
	}
	printf("\n");
//prop0.erase(prop0.begin() + 4);
//prop0.erase(prop0.begin() + 2);
//prop0.erase(prop0.begin() + 0);
//nnn.erase(nnn.begin() + 4);
//nnn.erase(nnn.begin() + 2);
//nnn.erase(nnn.begin() + 0);
//	_propErrorCorrect(prop0, nnn, response, glabel, 1);
	printf("pN=%d(%d): ", (int)nnn.size(), (int)prop0.size());for(i=0;i<(int)nnn.size();i++) printf("%d,",(int)nnn[i]);printf("\n");
	pN = (int)prop0.size();
	for(i = 0; i < pN; i++)
	{	printf("%d: ", i);
		for(j = 0; j < mylik.clustersz; j++) printf("%d,",(int)prop0[i][j]); printf("\n");
	}
	printf("\n");

/*	vector<int> P(mylik.clustersz), R, X;
	vector<vector<int> > clique;
	for(i = 0; i < mylik.clustersz; i++) P[i] = i;
	BronKerbosch1(score, P, R, X, clique, -6);

	vector<vector<int> > gcluster;
	vector<int> gmap(mylik.clustersz, -1);
	while((int)clique.size() > 0)
	{	vector<MySortType> cliqueW((int)clique.size());
		for(i = 0; i < (int)cliqueW.size(); i++)
		{	cliqueW[i].index = i;
			cliqueW[i].score = (double)clique[i].size();
			for(j = 0; j < (int)clique[i].size() - 1; j++)
				for(k = j + 1; k < (int)clique[i].size(); k++)
					cliqueW[i].score *= 1. - 1. / (1. + exp(score[clique[i][j]][clique[i][k]]));
		}
		sort(cliqueW.begin(), cliqueW.end());
		k = cliqueW[(int)cliqueW.size() - 1].index;
		vector<int> t = clique[k];
		for(i = 0; i < (int)t.size(); i++) gmap[t[i]] = (int)gcluster.size();
		gcluster.push_back(t);
		clique.erase(clique.begin() + k);
		for(i = (int)clique.size() - 1; i >= 0; i--)
		{	for(j = (int)clique[i].size() - 1; j >= 0; j--)
			{	for(k = 0; k < (int)t.size(); k++)
					if(t[k] == clique[i][j]) break;
				if(k < (int)t.size()) clique[i].erase(clique[i].begin() + j);
			}
			if((int)clique[i].size() == 0) clique.erase(clique.begin() + i);
		}
	}
	printf("{\n");
	for(i = 0; i < (int)gcluster.size(); i++)	
	{	for(j = 0; j < (int)gcluster[i].size(); j++)
			printf("%d,",gcluster[i][j]);
		printf("\n");
	}
	printf("}\n");

	vector<int> lab(pN, -1);
	for(i = 0; i < pN; i++)
	{	for(j = 0; j < mylik.clustersz; j++)
			if(prop0[i][j] > 0) 
			{	if(lab[i] < 0) lab[i] = gmap[j];
				else if(lab[i] != gmap[j]) { lab[i] = -1; break; }
			}
	}
	vector<double> ttnn;
	vector<vector<double> > ttpp;
	vector<int> pmap(mylik.clustersz, -1);
	for(i = 0; i < pN; i++)
	{	if(lab[i] < 0)
		{	ttnn.push_back(nnn[i]);
			ttpp.push_back(prop0[i]);
		}
		else 
		{	if(pmap[lab[i]] < 0)
			{	pmap[lab[i]] = (int)ttpp.size();
				ttnn.push_back(0);
				ttpp.push_back(vector<double>(mylik.clustersz, 0));
			}
			j = pmap[lab[i]];
			for(k = 0; k < mylik.clustersz; k++) ttpp[j][k] = (ttpp[j][k] * ttnn[j] + prop0[i][k] * nnn[i]) / (ttnn[j] + nnn[i]);
			ttnn[j] += nnn[i];
		}
	}
	prop0 = ttpp;
	nnn = ttnn;
	printf("pN=%d-->%d ", pN, (int)prop0.size());
	pN = (int)prop0.size();
*/

vector<double> tn = nnn;
sort(tn.begin(), tn.end());
for(i = 0; i < (int)tn.size(); i++) printf("%d,", (int)tn[i]);
printf("\n");
	vector<vector<double> > newprop;
	vector<double> newnnn;
	_groupComposition(prop0, nnn, newprop, newnnn, 1., 1.);
	printf("pN=%d --> %d\n", pN, (int)newprop.size());
tn = newnnn;
sort(tn.begin(), tn.end());
for(i = 0; i < (int)tn.size(); i++) printf("%d,", (int)tn[i]);
printf("\n");
	prop0 = newprop;
	nnn = newnnn;
	pN = (int)prop0.size();
	for(i = 0; i < pN; i++)
	{	printf("%d: ", i);
		for(j = 0; j < mylik.clustersz; j++) printf("%d,",(int)prop0[i][j]); printf("\n");
	}
	printf("\n");

if(false)
{
vector<vector<double> > newprop;
vector<double> newnnn;
for(i=0;i<pN;i++)for(j=0;j<mylik.clustersz;j++) prop0[i][j]=exp(prop0[i][j]);
_groupComposition1(prop0, nnn, newprop, newnnn);
printf("pN=%d-->",pN);
prop0 = newprop;
nnn = newnnn;
pN = (int)prop0.size();
for(i=0;i<pN;i++)for(j=0;j<mylik.clustersz;j++) prop0[i][j]=log(prop0[i][j]+1e-100);
printf("%d\n", pN);
}

FILE *fff = fopen("t.txt","w");
for(i = 0; i < pN; i++)
{	fprintf(fff,"%f ", i, nnn[i]);
	for(j = 0; j < mylik.clustersz; j++) fprintf(fff,"%f ", prop0[i][j]);
	fprintf(fff,"\n");
}
fclose(fff);

	vector<vector<double> > cprop0 = prop0;
	for(i = 0; i < pN; i++)
	{	for(j = 0; j < mylik.clustersz; j++)
		{	cprop0[i][j] = exp(prop0[i][j]);
			if(j > 0) cprop0[i][j] += cprop0[i][j - 1];
		}
	}
	vector<double> ncumsum(pN, 0);
	for(i = 0; i < pN; i++) 
	{	ncumsum[i] = nnn[i] / totalcn;
		nnn[i] = log(ncumsum[i]);
		ncumsum[i] = pow((exp(nnn[i])*totalcn),0.9);
		if(i > 0) ncumsum[i] += ncumsum[i - 1];
	}


/*
{
FILE *ff = fopen("t.txt", "w");
for(i = 0; i < (int)prop0.size(); i++)
{	for(j = 0; j < (int)prop0[i].size(); j++)
		fprintf(ff, "%f ", prop0[i][j]);
	fprintf(ff, "%d\n", (int)(exp(nnn[i]) * (double)totalcn));
}
fclose(ff);
}
*/
//int ccccc = 0;

	vector<MySortType> pp(mydata.L), ppt(mydata.L * 10);
        for(i = 0; i < (int)ppt.size(); i++)
        {       
//double uuuuu = gsl_ran_flat(gammar, 0, ncumsum[pN - 1]);
//for(k = 0; k < pN; k++) if(uuuuu <= ncumsum[k]) break;
//int kkkkk = k;
		int kkk = -1;
		bool renew = true;
                vector<double> lp0(pN, 0), tlp1 = lp0, wwlp0 = lp0, wwlp1 = lp0, twlp1 = lp0;
		double wlp1 = 0, wlp0 = 0, dwlp1 = 0;
                double lp1 = 0, maxlp, s, un;
                float yp[maxymsz * mydata.totalN], xp[maxxmsz * mydata.totalN + 1], *ypp, *xpp;
		ypp = &yp[0]; xpp = &xp[0];
		for(int jj = 0; jj < mydata.indN; jj++)
                {       j = groupid[jj].index;
			vector<double> rt;
			int l, ll, lll, idst = mydata.indIndex[j], ided = mydata.indIndex[j + 1];
			if(i < mydata.L)
			{	lll = 0;
				for(k = idst; k < ided; k++)
				{	ll = i * ymsz[k];
					for(l = 0; l < ymsz[k]; l++)
						ypp[lll++] = dataYP[k][ll++];
				}
				if(maxxmsz > 0)
				{	lll = 0;
					for(k = idst; k < ided; k++)
					{	ll = i * xmsz[k];
						for(l = 0; l < xmsz[k]; l++)
						{	xpp[lll++] = dataXP[k][ll++];
						}
					}
				}
				_getLpVar(ypp, xpp, idst, ided, ymsz + idst, xmsz + idst, prop0, tpara.priorW, mylik.clustersz, rt);
        	                for(k = 0; k < pN; k++) 
				{	lp0[k] += rt[k];
					tlp1[k] += rt[k];
				}
				if(jj ==  mydata.indN - 1 || groupid[jj].score != groupid[jj+1].score)
				{	maxlp = tlp1[0] + nnn[0]; s = 0;
					for(k = 1; k < pN; k++) maxlp = max(maxlp, tlp1[k] + nnn[k]);
					for(k = 0; k < pN; k++) 
					{	s += exp(tlp1[k] + nnn[k] - maxlp); 
						tlp1[k] = 0;
					}
					lp1 += log(s) + maxlp;
				}
			}
			
			if(jj==0 || groupid[jj].score != groupid[jj-1].score)
			{	if(renew)
				{	un = gsl_ran_flat(gammar, 0, ncumsum[pN - 1]);
					for(kkk = 0; kkk < pN; kkk++) if(un <= ncumsum[kkk]) break;
				}
				un = gsl_ran_flat(gammar, 0., 1.);
				if(un <= 1. - pow(1. - 0.5, 1. / (double)rN)) renew = true;
				else renew = false;
			}
			k = kkk;
//k=kkkkk;
			int ok = k;
			int sel[ided - idst], selk = k;
			int llly = 0, lllx = 0;
			for(ll = 0; ll < ided - idst; ll++)
			{	un = gsl_ran_flat(gammar, 0, cprop0[k][mylik.clustersz - 1]);	
				for(l = 0; l < mylik.clustersz; l++) if(un <= cprop0[k][l]) break;
				mylik.simData(1, idst + ll, l, ypp + llly);
				sel[ll] = l;
				for(k = 0; k < xmsz[idst + ll]; k++) xpp[lllx + k] = 0;
				llly += ymsz[idst + ll];
				lllx += xmsz[idst + ll];
			}
			_getLpVar(ypp, xpp, idst, ided, ymsz + idst, xmsz + idst, prop0, tpara.priorW, mylik.clustersz, rt);
                        for(k = 0; k < pN; k++) 
			{	wwlp0[k] += rt[k];
				wwlp1[k] += rt[k];
				twlp1[k] += rt[k];
			}
                	
			if(jj == mydata.indN - 1 || groupid[jj].score != groupid[jj+1].score)
			{
				double tmaxlp = twlp1[0] + log(ncumsum[0] / ncumsum[pN-1]), ts = 0; 
				maxlp = wwlp1[0] + nnn[0]; s = 0;
				for(k = 1; k < pN; k++)
				{	tmaxlp = max(tmaxlp, twlp1[k] + log((ncumsum[k]-ncumsum[k-1])/ncumsum[pN-1]));
					maxlp = max(maxlp, wwlp1[k] + nnn[k]);
				}
				for(k = 0; k < pN; k++) 
				{	if(k==0) ts += exp(twlp1[k] + log(ncumsum[k]/ncumsum[pN-1]) - tmaxlp);
					else ts += exp(twlp1[k] + log((ncumsum[k]-ncumsum[k-1])/ncumsum[pN-1]) - tmaxlp);
					if(jj == mydata.indN - 1 || renew) twlp1[k] = 0;
					s += exp(wwlp1[k] + nnn[k] - maxlp);
					wwlp1[k] = 0;
				}
				wlp1 += log(s) + maxlp;
				if(jj == mydata.indN - 1 || renew) dwlp1 += log(ts) + tmaxlp;
			}

			ypp += llly;
			xpp += lllx;
		}       
		double maxlp0, s0, tmax;
		if(i < mydata.L)
 		{	maxlp0 = lp0[0] + nnn[0] * (double)rN; s0 = 0;
			for(j = 1; j < pN; j++) maxlp0 = max(maxlp0, lp0[j] + nnn[j] * (double)rN);
			for(j = 0; j < pN; j++) s0 += exp(lp0[j] + nnn[j] * (double)rN - maxlp0);
			tmax = max(lp1, maxlp0);
			lp1 = log(max(exp(lp1 - tmax) - s0 * exp(maxlp0 - tmax),1e-100)) + tmax;

			maxlp = lp0[0] + nnn[0]; s = 0;
			for(j = 1; j < pN; j++) maxlp = max(maxlp, lp0[j] + nnn[j]);
			for(j = 0; j < pN; j++) s += exp(lp0[j] + nnn[j] - maxlp);
			s += s0 * exp(maxlp0-maxlp);
			s = log(s) + maxlp;
			maxlp = max(maxlp, lp1);maxlp=max(maxlp,s);
			vprob[i] = lp1 - s;//exp(lp1 - maxlp) / (exp(lp1 - maxlp) + exp(s - maxlp));
			pp[i].index = i;
			pp[i].score = lp1-s;
		}

 		maxlp0 = wwlp0[0] + nnn[0] * (double)rN; s0 = 0;
		for(j = 1; j < pN; j++) maxlp0 = max(maxlp0, wwlp0[j] + nnn[j] * (double)rN);
		for(j = 0; j < pN; j++) s0 += exp(wwlp0[j] + nnn[j] * (double)rN - maxlp0);
		tmax = max(wlp1, maxlp0);
		wlp1 = log(max(exp(wlp1 - tmax) - s0 * exp(maxlp0 - tmax),1e-100)) + tmax;

                maxlp = wwlp0[0] + nnn[0]; s = 0;
		for(j = 1; j < pN; j++) maxlp = max(maxlp, wwlp0[j] + nnn[j]);
		for(j = 0; j < pN; j++) s += exp(wwlp0[j] + nnn[j] - maxlp);
		wlp0 = log(s) + maxlp;
		s += s0 * exp(maxlp0-maxlp);
		s = log(s) + maxlp;
		maxlp = max(maxlp, wlp1);maxlp = max(maxlp, s);

		ppt[i].index = i;
		ppt[i].score = //fabs(yp[0]-yp[1]);
				wlp1-s;
				//exp(wlp1 - maxlp) / (exp(wlp1 - maxlp) + exp(s - maxlp));
		ppt[i].weight = exp(wlp0 - dwlp1);//k;
//printf("%d %f %f | %f %f %f | %f\n", i, yp[0], yp[1], wlp0, dwlp1, ppt[i].weight, ppt[i].score);
//printf("c %f %d\n",vprob[i], i);fflush(stdout);

//if(ppt[i].score >= 0.9) ccccc++;
//ppt[i].weight = 1.;

//pN -= rN;
	}

	sort(pp.begin(), pp.end());
	sort(ppt.begin(), ppt.end());
	reverse(pp.begin(), pp.end());	
	reverse(ppt.begin(), ppt.end());
	double ss = 0, ff, off = 1.;
	j = 0;
	for(i = 0; i < mydata.L; i++)
	{	while(j < (int)ppt.size() && ppt[j].score >= pp[i].score) { ss += ppt[j].weight; j++; }
		vprob[pp[i].index + mydata.L] = min(1., ss / (double)ppt.size());
	}
	for(i = mydata.L - 1; i >= 0; i--)
	{	ff = vprob[pp[i].index + mydata.L] * (double)(mydata.L) / (double)(i + 1);
		off = vprob[pp[i].index + mydata.L*2] = min(off, ff);
	}
	
}

void genomicTensor::_groupComposition(vector<vector<double> > &ggg, vector<double> &nnn, vector<vector<double> > &prop, vector<double> &newnnn, double a, double b)
{	int i, j, k, l, N = (int)ggg.size(), L = (int)ggg[0].size();

	vector<vector<double> > score;
	_pairEnrich(ggg, nnn, score);

	double buff  = 3.;
	vector<vector<double> > sim(N, vector<double>(N, -100000000.));
	for(i = 0; i < N; i++)
	{	for(j = i; j < N; j++)
		{	double f = 0, n = 0;
			for(k = 0; k < L; k++)
				if(ggg[i][k] > 0)
				{	for(l = 0; l < L; l++)
						if(ggg[j][l] > 0)
						{	f += min(0.,score[k][l]);//log(1. - 1. / (1. + exp(score[k][l])));
							n++;
						}
				}
			sim[i][j] = sim[j][i] = f + buff;
		}
	}
for(i=0;i<N;i++)
{	for(j=0;j<N;j++)printf("%3.3f ",sim[i][j]);printf("\n");
}
/*
for(i = N - 1; i >= 0; i--)
	if(sim[i][i] < 0)
	{	ggg.erase(ggg.begin() + i);
		nnn.erase(nnn.begin() + i);
		sim.erase(sim.begin() + i);
		for(j = 0; j < (int)sim.size(); j++) sim[j].erase(sim[j].begin() + i);
	}
printf("N=%d-->%d\n", N, (int)ggg.size());
N = (int)ggg.size();
*/
a=1e-10;b=1e-10;
	double A = a * (double)L, B = b * N;
	double *GGG = new double[N * L * 2 + N * 2], *NNN = GGG + N * L, *PROP = NNN + N, *NEWNNN = PROP + N * L;
	double *tggg, *tprop;
	double totaln = 0;
	tggg = GGG; tprop = PROP;
	for(i = 0; i < N; i++)
	{	for(j = 0; j < L; j++, tggg++, tprop++)
		{	(*tggg) = ggg[i][j];
			NNN[i] = nnn[i];
			(*tprop) = (*tggg) * NNN[i]; 
		}
		totaln += NNN[i];
	}
	int gN = 0;
	int cn[N], map[N], nn[N];
	newnnn.resize(N);
	tggg = GGG;
	for(i = 0; i < N; i++)
	{	k = 0; l = 0;
		for(j = 0; j < L; j++, tggg++) k += (int)(*tggg);//, l+=(int)prop[i][j];
		NEWNNN[i] = NNN[i] * k;
		cn[i] = (int)NEWNNN[i];
		nn[i] = (int)NNN[i];
		map[i] = i;
		if(nnn[i] > 0) gN++;
	}
	
	double *gLP = new double[N * N * 2], *TP1 = gLP + N * N, *gp, *gm;
	double sumLP[N], F0[N];
	tggg = GGG;
	gm = gLP;
	for(l = 0; l < N; l++, tggg += L)
	{	double sum = 0, f0 = 0;
		for(k = 0; k < L; k++)
		{	f0 -= gsl_sf_lngamma(tggg[k]+1);
			sum += tggg[k];
		}
		f0 += gsl_sf_lngamma(sum+1);
		F0[l] = f0;
		for(k = 0; k < N; k++, gm++) *gm = f0;
		for(k = 0; k < L; k++)
		{	gp = gLP + l * N;
			tprop = PROP + k;
			for(j = 0; j < N; j++, gp++, tprop += L)
				(*gp) += log((*tprop + a) / (cn[j] + A)) * tggg[k];
		}
	}
	gp = gLP;
	for(k = 0; k < N; k++, gp += N)
	{	double s = 0, maxlp = *gp + log((nn[0] + b) / (totaln + B));
		for(l = 1; l < N; l++) maxlp = max(maxlp, gp[l] + log((nn[l] + b) / (totaln + B)));
		for(l = 0; l < N; l++) s += exp(gp[l] + log((nn[l] + b) / (totaln + B)) - maxlp);
		sumLP[k] = log(s) + maxlp;
	}

	int burnin = 0, mcmc = 50;
	for(i = 0; i < burnin + mcmc; i++)
	{	bool flag = false;
		printf("%d ", i);fflush(stdout);
		for(j = 0; j < N; j++) printf("%d:%d ", j,nn[j]);printf("\n");
		tggg = GGG;
		for(j = 0; j < N; j++, tggg += L)
		{	printf("%d:%d ", j,map[j]);fflush(stdout);
			int mj = map[j];
			tprop = PROP + mj * L;
			for(k = 0; k < L; k++, tprop++)
				(*tprop) -= tggg[k] * NNN[j];
			cn[mj] -= (int)NEWNNN[j];
			nn[mj] -= (int)NNN[j];
			double lp[N], lp0 = MINUSINFINITE;
			double maxlp = (double)MINUSINFINITE;
			
			double *tm = TP1 + mj, *tp, *ttgg = GGG;
			for(k = 0; k < N; k++, tm += N)
			{	*tm = F0[k];
				tprop = PROP + mj * L;
				for(l = 0; l < L; l++, ttgg++, tprop++)
					if(*ttgg > 0) (*tm) += log((*tprop + a) / (cn[mj] + A)) * (*ttgg);
			}
			for(k = 0; k < N; k++)
			{	lp[k] = (double)MINUSINFINITE;
				if(cn[k] > 0 || (cn[k] == 0 && lp0 <= MINUSINFINITE))
				{	lp[k] = 0;
			/*		for(l = 0; l < L; l++)
					{	if(GGG[j*N+l] > 0.01)
							lp[k] += gsl_sf_lngamma(PROP[k*N+l] + GGG[j*N+l] * NNN[j] + a) - gsl_sf_lngamma(PROP[k*N+l] + a);
					}
					lp[k] -= gsl_sf_lngamma((double)cn[k] + NEWNNN[j] + A) - gsl_sf_lngamma((double)cn[k] + A); 
			*/		
					if(k != mj)
					{	tp = TP1 + k; ttgg = GGG;
						for(int kk = 0; kk < N; kk++, tp += N)
						{	*tp = F0[kk];
							tprop = PROP + k * L;
							for(l = 0; l < L; l++, ttgg++, tprop++)
								if(*ttgg > 0) (*tp) += log((*tprop + tggg[l] * NNN[j] + a) / (cn[k] + NEWNNN[j] + A)) * (*ttgg);
						}
					}
					gp = gLP + k; tp = TP1 + k;
					gm = gLP + mj; tm = TP1 + mj;
					for(l = 0; l < N; l++, gp += N, gm += N, tp += N, tm += N)
					{	double newsumLP = sumLP[l], tmax = max(newsumLP, *tp);
						if(k != mj)
						{	newsumLP = log(max(exp(newsumLP - tmax) - exp(*gp + log((nn[k] + b) / (totaln + B)) - tmax) + exp(*tp + log((nn[k] + NNN[j] + b) / (totaln + B)) - tmax), 1e-100)) + tmax;
							tmax = max(newsumLP, *tm);
							newsumLP = log(max(exp(newsumLP - tmax) - exp(*gm + log((nn[mj] + NNN[j] + b) / (totaln + B)) - tmax) + exp(*tm + log((nn[mj] + b) / (totaln + B)) - tmax), 1e-100)) + tmax;
						}
						lp[k] += newsumLP * NNN[l];
					}
			
//printf("%d(%d).%d:(%d) %f->",j,map[j],k,cn[k],lp[k]);double dd = lp[k];
				double bb=totaln / 1.;	
					lp[k] += gsl_sf_lngamma((double)nn[k] + nnn[j] + bb) - gsl_sf_lngamma((double)nn[k] + bb); 
			//	int remove = (int)(nn[mj] == 0), add = (int)(nn[k] == 0);
			//	lp[k]+=(double)(remove-add)*10.*(double)(gN+add-remove > 10);
//printf("%f (%f)| ",lp[k],lp[k]-dd);
				for(int tt = 0; tt < N; tt++) 
				{	//if(tt != j && map[tt] == k && sim[tt][j] < 0) { lp[k] += sim[tt][j] * 10.; }
				//	if((map[tt] == k || tt == j) && sim[tt][tt] < 0 && nn[k] <= nnn[tt]) lp[k] += sim[tt][tt];//1000000; 
				}
				//if(sim[j][j] < 0 && nn[k]  < nnn[j]) lp[k] -= 1000000.;

					maxlp = max(maxlp, lp[k]);
					if(cn[k] == 0) lp0 = lp[k];
				}
			}
			if(i < burnin)
			{	for(k = 0; k < N; k++)
				{	lp[k] = exp(lp[k] - maxlp);
				}
			}
			k = _sample(lp, N, (int)(i>=burnin));
//printf("%d\n",k);fflush(stdout);
			if(mj != k && i >= burnin) //printf("%d %d:%d->%d %f %f, %d %d\n",i,j,map[j],k, lp[map[j]], lp[k], nn[map[j]], nn[k]),
				flag = true;
			if(mj != k)
			{	gp = gLP + k; tp = TP1 + k;
				gm = gLP + mj; tm = TP1 + mj;
				for(l = 0; l < N; l++, gp += N, gm += N, tp += N, tm += N)
				{	double newsumLP = sumLP[l], tmax = max(newsumLP, *tp);
					newsumLP = log(max(exp(newsumLP - tmax) - exp(*gp + log((nn[k] + b) / (totaln + B)) - tmax) + exp(*tp + log((nn[k] + NNN[j] + b) / (totaln + B)) - tmax), 1e-100)) + tmax;
double os = newsumLP;
					tmax = max(newsumLP, *tm);
					newsumLP = log(max(exp(newsumLP - tmax) - exp(*gm + log((nn[mj] + NNN[j] + b) / (totaln + B)) - tmax) + exp(*tm + log((nn[mj] + b) / (totaln + B)) - tmax), 1e-100)) + tmax;
					sumLP[l] = newsumLP;
				}
				gp = gLP + k; tp = TP1 + k;
				gm = gLP + mj; tm = TP1 + mj;
				for(l = 0; l < N; l++, gp += N, gm += N, tp += N, tm += N) 
				{	*gp = *tp;
					*gm = *tm;
				}
			}
			if(nn[mj]==0) gN--; if(nn[k]==0) gN++;
			mj = map[j] = k;
			tprop = PROP + mj * L;
			for(k = 0; k < L; k++, tprop++)
				*tprop += tggg[k] * NNN[j];
			cn[mj] += (int)NEWNNN[j];
			nn[mj] += (int)NNN[j];
		}
		if(!flag) break;
	}
	delete gLP;
	
	vector<double> totalprop(L,10);
	tprop = PROP;
	for(i = 0; i < N; i++)
		for(j = 0; j < L; j++, tprop++)
			totalprop[j] += (*tprop + a);
	double totalcn = 0;
	for(i = 0; i < L; i++) totalcn += totalprop[i];
	
	tprop = PROP;
	for(i = 0; i < N; i++, tprop += L)
	{	for(j = 0; j < L; j++) tprop[j] = log((tprop[j] + a) / (cn[i] + A));
						//log((tprop[j] + a * totalprop[j]/totalcn)/(cn[i]+A));
		NEWNNN[i] = (double)nn[i] + b;
	}

	prop.clear();
	newnnn.clear();
	tprop = PROP;
	for(i = 0; i < N; i++, tprop += L)
	{	if(NEWNNN[i] > b + 10)
		{	vector<double> row(L);
			for(j = 0; j < L; j++) row[j] = tprop[j];
			prop.push_back(row);
			newnnn.push_back(NEWNNN[i]);
		}
	}
	
	delete GGG;
}

void genomicTensor::_pairEnrich(vector<vector<double> > const &prop, vector<double> const &nnn, vector<vector<double> > &score)
{	int i, j, k;
	int K = (int)prop[0].size();
	double sum = 0;
	vector<vector<double> > A(K, vector<double>(K, 0));
	for(i = 0; i < (int)prop.size(); i++)
	{	sum = 1e-100;
		for(j = 0; j < K; j++) sum += prop[i][j];
		for(j = 0; j < K; j++)
			for(k = j; k < K; k++)
				if(prop[i][j] > 0 && prop[i][k] > 0)
				{	double t = prop[i][j] * (prop[i][k] - (double)(j==k)) / (1. + (double)(j==k)) * nnn[i] / (sum - 1.);
					A[j][k] += t;
					A[k][j] += t;
				}
	}
	sum = 1e-100;
	vector<double> P(K, 0);
	for(i = 0; i < K; i++)
	{	for(j = 0; j < K; j++) P[i] += A[i][j];
		sum += P[i];
	}
	for(i = 0; i < K; i++) P[i] /= sum;
	vector<vector<double> > E(K, vector<double>(K));
	for(i = 0; i < K; i++)
	{	for(j = i; j < K; j++)
			E[i][j] = E[j][i] = P[i] * P[j] * sum;
	}
	score = A;
	for(i = 0; i < K; i++)
	{	for(j = i; j < K; j++)
		{	score[i][j] = score[j][i] = (A[i][j] - E[i][j]) / sqrt(E[i][j] + 1.);
		}
	}
}

void genomicTensor::BronKerbosch1(vector<vector<double> > const &score, vector<int> P, vector<int> R, vector<int> X, vector<vector<int> > &clique, double scorecut)
{	int i, j, k, l;

/*			printf("[P: ");
			for(i = 0; i < (int)P.size(); i++) printf("%d,",P[i]);
			printf("]--");
			printf("[R: ");
			for(i = 0; i < (int)R.size(); i++) printf("%d,",R[i]);
			printf("]--");
			printf("[X: ");
			for(i = 0; i < (int)X.size(); i++) printf("%d,",X[i]);
			printf("]\n");
*/
	if((int)P.size() == 0 && (int)X.size() == 0) 
	{	if((int)R.size() > 0)
		{	clique.push_back(R);
			printf("[ ");
			for(i = 0; i < (int)R.size(); i++) printf("%d,",R[i]);
			printf("]\n");
		}
		return;
	}
	
	while((int)P.size() > 0)
	{	int v = P[0];
		vector<int> nP, nX, nR = R;
		nR.push_back(v);
		for(i = 1; i < (int)P.size(); i++)
		{	if(score[v][P[i]] > scorecut) nP.push_back(P[i]);
		}	
		for(i = 0; i < (int)X.size(); i++)
		{	if(score[v][X[i]] > scorecut) nX.push_back(X[i]);
		}	
		BronKerbosch1(score, nP, nR, nX, clique, scorecut);
		P.erase(P.begin());
		X.push_back(v);
	}
}

void genomicTensor::_groupComposition1(vector<vector<double> > const &ggg, vector<double> const &nnn, vector<vector<double> > &prop, vector<double> &newnnn)
{
	vector<vector<double> > score;
	_pairEnrich(ggg, nnn, score);

	int N = (int)ggg.size(), L = (int)ggg[0].size();
	int i, j, k, l;
	double buff  = 3.;
	vector<vector<double> > sim(N, vector<double>(N, -100000000.));
	for(i = 0; i < N - 1; i++)
	{	for(j = i; j < N; j++)
		{	double f = 0, n = 0;
			for(k = 0; k < L; k++)
				if(ggg[i][k] > 0)
				{	for(l = 0; l < L; l++)
						if(ggg[j][l] > 0)
						{	f += min(0., score[k][l]);//log(1. - 1. / (1. + exp(score[k][l])));
							n++;
						}
				}
			sim[i][j] = sim[j][i] = f + buff;
		}
	}
	
	vector<int> list(N);
	prop = ggg;
	for(i = 0; i < N; i++)
	{	for(j = 0; j < L; j++) prop[i][j] *= nnn[i];
		list[i] = i;
	}
	newnnn = nnn;
	bool a_inner = true;
	do
	{	int tN = (int)list.size();
		vector<MySortType> me(tN * (tN - 1) / 2);
		k = 0;
		for(i = 0; i < tN - 1; i++)
			for(j = i + 1; j < tN; j++)
			{	me[k].score = sim[list[i]][list[j]];
				me[k].index = i * tN + j;
				k++;
			}
		sort(me.begin(), me.end());
		if(me[(int)me.size() - 1].score > 0)
		{	i = me[(int)me.size() - 1].index;
			j = i % tN;
			i = (int)(i / tN);
			for(k = 0; k < L; k++) 
			{	prop[list[i]][k] += prop[list[j]][k];
				prop[list[j]][k] = 0;
			}
			newnnn[list[i]] += newnnn[list[j]];
			newnnn[list[j]] = 0;
			list.erase(list.begin() + j);
			int ii = list[i];
			for(k = 0; k < tN - 1; k++)
				if(k != i)
				{	int jj = list[k];
					double f = 0, n = 0;
					for(int kk = 0; kk < L; kk++)
						if(prop[ii][kk] > 0)
						{	for(l = 0; l < L; l++)
								if(prop[jj][l] > 0)
								{	f += log(1. - 1. / (1. + exp(score[kk][l])));
									n++;
								}
						}	
					sim[ii][jj] = sim[jj][ii] = f + buff;	
				}
		}
		else a_inner = false;
	} while(a_inner);

	j = 0;
	for(i = 0; i < N; i++)
		if(newnnn[i] > 0)
		{	prop[j] = prop[i];
			for(k = 0; k < L; k++) prop[j][k] = ((prop[i][k] + 0) / newnnn[i]);
			newnnn[j] = newnnn[i];
			j++;
		}
	prop.resize(j);
	newnnn.resize(j);
}

void genomicTensor::_propErrorCorrect(vector<vector<double> > &ggg, vector<double> &nnn, vector<int> const &response, vector<vector<int> > &glabel, double a)
{	int i, j, k, l, pN = (int)ggg.size();
	vector<vector<double> > prop0 = ggg;
	for(i = 0; i < pN; i++)
	{	double s = 0;
		for(j = 0; j < mylik.clustersz; j++) s += ggg[i][j] + 1e-100;
		for(j = 0; j < mylik.clustersz; j++) prop0[i][j] = log((ggg[i][j] + 1e-100)/s);
	}
	vector<MySortType> groupid(mydata.indN);
	for(i = 0; i < mydata.indN; i++)
	{	groupid[i].index = i;
		groupid[i].score = (double)response[i];
	}
	sort(groupid.begin(), groupid.end());
	int rN = 0;
	for(i = 0; i < mydata.indN; i++)
		if(i == mydata.indN - 1 || groupid[i].score != groupid[i+1].score) rN++;

	vector<vector<double> > score;
	_pairEnrich(ggg, nnn, score);

	double buff  = 6.;
	vector<vector<double> > sim(pN, vector<double>(pN, -100000000.));
	for(i = 0; i < pN; i++)
	{	for(j = i; j < pN; j++)
		{	double f = 0, n = 0;
			for(k = 0; k < mylik.clustersz; k++)
				if(ggg[i][k] > 0)
				{	for(l = 0; l < mylik.clustersz; l++)
						if(ggg[j][l] > 0)
						{	f += min(0.,score[k][l]);//log(1. - 1. / (1. + exp(score[k][l])));
							n++;
						}
				}
			sim[i][j] = sim[j][i] = f + buff;
		}
	}
for(i=0;i<pN;i++)printf("%f ",sim[i][i]);

        float yp[maxymsz * mydata.totalN], xp[maxxmsz * mydata.totalN + 1];
	int burnin = 0, mcmc = 50;
	for(int r = 0; r < burnin + mcmc; r++)
	{	bool flag = false;
		for(i = 0; i < mydata.L; i++)
		{	vector<double> lp(pN, 0);
			int rid = 0;
			for(int jj = 0; jj < mydata.indN; jj++)
			{	j = groupid[jj].index;
				vector<double> rt;
				int ll, lll, idst = mydata.indIndex[j], ided = mydata.indIndex[j + 1];
				lll = 0;
				for(k = idst; k < ided; k++)
				{	ll = i * ymsz[k];
					for(l = 0; l < ymsz[k]; l++)
						yp[lll++] = dataYP[k][ll++];
				}
				if(maxxmsz > 0)
				{	lll = 0;
					for(k = idst; k < ided; k++)
					{	ll = i * xmsz[k];
						for(l = 0; l < xmsz[k]; l++)
						{	xp[lll++] = dataXP[k][ll++];
						}
					}
				}
				_getLpVar(yp, xp, idst, ided, ymsz + idst, xmsz + idst, prop0, tpara.priorW, mylik.clustersz, rt);
        	                for(k = 0; k < pN; k++) 
				{	lp[k] += rt[k];
				}
				if(jj ==  mydata.indN - 1 || groupid[jj].score != groupid[jj+1].score) 
				{	double tlp[pN], maxlp;
if(i>=(int)glabel.size() || rid >= (int)glabel[0].size() || glabel[i][rid] < 0 || glabel[i][rid] >= pN)
{	printf("!!!%d,%d,%d,%d,%d,%d\n",i,(int)glabel.size(), rid,(int)glabel[0].size(),glabel[i][rid],pN);fflush(stdout);
}
					for(k = 0; k < pN; k++)
					{	tlp[k] = lp[k] + log(nnn[k] - (double)(k == glabel[i][rid]) + a);
						if(k == 0 || maxlp < lp[k]) maxlp = tlp[k];
					}
					if(r < burnin)
					{	for(k = 0; k < pN; k++) tlp[k] = exp(tlp[k] - maxlp);	
					}
					k = _sample(tlp, pN, (int)(r >= burnin));
int t=glabel[i][rid];if(sim[t][t]>=0) k=glabel[i][rid];
					if(k != glabel[i][rid])
					{	nnn[glabel[i][rid]]--;
						glabel[i][rid] = k;
						nnn[k]++;
						flag = true;
					}
					lp.clear();
					lp.resize(pN, 0);
					rid++;
				}
			}
		}
		if(!flag && r >= burnin) break;
	}
/*
	vector<vector<double> > newprop(pN, vector<double>(mylik.clustersz, 0));
	for(i = 0; i < mydata.L; i++)
	{	int rid = 0;
		for(int jj = 0; jj < mydata.indN; jj++)
		{	j = groupid[jj].index;
			k = glabel[i][rid];
			int idst = mydata.indIndex[j], ided = mydata.indIndex[j+1];
			for(l = idst; l < ided; l++)
				newprop[k][(int)mydata.data[j * mydata.L + i]]++;
			if(jj ==  mydata.indN - 1 || groupid[jj].score != groupid[jj+1].score) rid++;
		}
	}
	for(i = 0; i < pN; i++)
	{	printf("[%d: ",i);
		for(j = 0; j < mylik.clustersz; j++)
		{	printf("%d,",(int)newprop[i][j]);newprop[i][j] /= (nnn[i] + 1e-100);
		}
		printf("] %d \n", (int)nnn[i]);
	}
	ggg = newprop;
*/
	for(i = (int)nnn.size() - 1; i >= 0; i--)
		if(nnn[i] == 0)
		{	ggg.erase(ggg.begin() + i);	
			nnn.erase(nnn.begin() + i);
		}
}	

void genomicTensor::_refreshParameter(bool updateparam)
{	int i, j;
	mylik.clearParameter();
	tpara.param.clear();
	vector<vector<float> > eee(tpara.maxK, vector<float>(mylik.clustersz, 0));
	tpara.param.resize(mydata.L, eee);
        float *dp = mydata.data;
	if(updateparam)
	{	for(int id = 0; id < mydata.indN; id++)
               		for(i = mydata.indIndex[id]; i < mydata.indIndex[id+1]; i++)
      	        	{	vector<vector<vector<float > > >::iterator ipm = tpara.param.begin();
				int *pp = tpara.pop + id * mydata.L;
                       		for(j = 0; j < mydata.L; j++, dp++, ipm++, pp++)
				{	_updatePara_add((*ipm)[(int)(*pp)], *dp, 1., i, j);
				}
               		}
	}
	else
	{	for(int id = 0; id < mydata.indN; id++)
               		for(i = mydata.indIndex[id]; i < mydata.indIndex[id+1]; i++)
      	        	{	vector<vector<vector<float > > >::iterator ipm = tpara.param.begin();
                       		for(j = 0; j < mydata.L; j++, dp++, ipm++)
				{	_updatePara_add((*ipm)[0], *dp, 1., i, j);
				}
               		}
	}
	mylik.updateLambda(tpara.priorW);
}

void genomicTensor::removeData(vector<vector<int> > const &list)
{
	int i, j, k;
	for(i = 0; i < (int)list.size(); i++)
		if((int)list[i].size() > 0)
		{	int keep[ymsz[i]];
			for(j = 0; j < ymsz[i]; j++) keep[j] = 0;
			printf("%d: ", i);
			for(j = 0; j < (int)list[i].size(); j++) { keep[list[i][j]] = 1; printf("%d,",list[i][j]); }; 
			printf("\n");
			k = 0; for(j = 0; j < ymsz[i]; j++) if(keep[j] == 0) keep[k++] = j;
			int omsz = ymsz[i];
			ymsz[i] = k;
			for(j = 0; j < k; j++) ymap[i][j] = ymap[i][keep[j]]; 
			float *dp = dataYP[i], *dp1 = dp; 
			for(j = 0; j < mydata.L; j++, dp1 += omsz)
			{	for(k = 0; k < ymsz[i]; k++)
				{	*dp = dp1[keep[k]];
					dp++;
				}
			}
			printf("%d:sz=%d ", i, ymsz[i]);
			for(j = 0; j < ymsz[i]; j++) printf("%d,", ymap[i][j]); printf("\n");
		}
//	for(i = 0; i < mydata.totalN; i++) printf("%f ", dataYP[i][3]); printf("\n");
//	exit(0);
//	for(i=0;i<mydata.totalN;i++) printf("%d:%d ", i, ymsz[i]);
}
