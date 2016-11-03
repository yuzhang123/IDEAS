#include "popTensor.h"

popTensor::popTensor(unsigned int rs):tensorHMMbase(rs)
{
	dataP = NULL;
}

popTensor::~popTensor()
{	if(dataP != NULL) delete dataP;
}

void popTensor::_getAsz(MYDATA const &mydata)
{	int i, j, maxK = 0;
	float const *p = mydata.data;
	for(i = 0; i < mydata.indN; i++)
		for(j = 0; j < mydata.L; j++, p++)
			maxK = max(maxK, (int)(*p) + 1); 
	for(i = 0; i < mydata.L; i++) mydata.asz[i] = maxK;
}

/////////////////////////////////////////////////////////////////////////////////////
void popTensor::run(char const *finput, int burnin, int mcmc, int maxHapK, int ploidity, double AA, char const *foutput, char const *fS, int admixN)
{
	int i;
	MYDATA mydata;
	
	trueploidity = mydata.ploidity = ploidity;
	vector<double> recomb;
	readInput(mydata, finput, fS, recomb);
	
	TENSORPARA tpara;
	initPara(tpara);
	tpara.rPrior = new double[mydata.L];
	for(i = 0; i < mydata.L; i++) 
		tpara.rPrior[i] = recomb[i] * (double)admixN / (double)max(admixN, 100);
	recomb.clear();

	tpara.priorW = 10.;
	tpara.maxHapK = maxHapK;
	tpara.burnin = burnin;
	tpara.mcmc = mcmc;
	if(mydata.indN * mydata.L > 10000000) tpara.imputedata = false;//memory cost can be too big, need to revise

	tpara.addsample = false;
	tpara.changeAlpha = -1;
	tpara.splitstop = 10;
	tpara.samplemaximum = true;
	if(AA > 0) tpara.probAlpha = AA;

	inferStructure(mydata, tpara, foutput, false);
	
	if(mydata.data != NULL) delete mydata.data;
	if(mydata.indIndex != NULL) delete mydata.indIndex;
}

/////////////////////////////////////////////////////////////////////////////////////
void popTensor::_updatePara_remove(vector<float> &para, float data, float b, int id, int pos)
{	para[(int)data] -= b;	
}

void popTensor::_updatePara_add(vector<float> &para, float data, float b, int id, int pos)
{	para[(int)data] += b;
}

//need to update pointers, bst and bed are breaks index, if they exist, otherwise they are positions
void popTensor::_sumforwardOne(int id, vector<int> const &breaks, int bst, int bed, TENSORPARA &tpara, int ploidity, float *&dp, float *&dp2, double const *&rp, int const *&nhp, vector<vector<vector<float> > >::const_iterator &ipm, int const *&aszp, double *&lpp, double *&lmm, double *&lss, int SSZ, int MK, vector<vector<int> > const &ss, vector<vector<int> > const &mlist, double *Vh2)
{	int i, j, k;

//printf("%d:%d(%d)-%d(%d) MK=%d ", id, breaks[bst], bst, breaks[bed], bed, MK);fflush(stdout);
	double tP1[MK * MK], tP2[MK * MK], mP1[MK], mP2[MK], sP1[MK], sP2[MK], cumsum[SSZ+4];
	int  path1[MK*MK*(bed-bst-2)/2], path2[MK*MK*(bed-bst-2)/2], stepn1, stepn2;

	stepn1 = _oneStrandProb(id, breaks, bst, bed, 0, tpara, rp, nhp, ipm, aszp, MK, &tP1[0], &sP1[0], &path1[0]);
	if(ploidity == 2) stepn2 = _oneStrandProb(id + 1, breaks, bst, bed, 1, tpara, rp, nhp, ipm, aszp, MK, &tP2[0], &sP2[0], &path2[0]);

/*
if(id==12 && bst==0)
{	printf("\n[%d:%d-%d], Vh=", id, breaks[bst], breaks[bed]);
	for(i = 0; i <= tpara.maxK; i++) printf("%f ", tpara.Vh[i]);
	printf("\ntP1=\n");
	for(i=0;i<MK;i++)
	{	for(j=0;j<MK;j++) printf("%f ", tP1[i*MK+j]);
		printf("\n");
	}
	printf("\ntP2=\n");
	for(i=0;i<MK;i++)
	{	for(j=0;j<MK;j++) printf("%f ", tP2[i*MK+j]);
		printf("\n");
	}
	printf("\n");
	printf("path1[%d]\n", stepn1);
	for(i=0;i<stepn1;i++)
	{	printf("[%d]\t", i);
		for(j=0;j<MK;j++)
		{	for(k=0;k<MK;k++)
				printf("%d ", path1[i*MK*MK+j*MK+k]);
			printf("\n\t");
		}
		printf("\n");
	}
	printf("path2[%d]\n", stepn2);
	for(i=0;i<stepn2;i++)
	{	printf("[%d]\t", i);
		for(j=0;j<MK;j++)
		{	for(k=0;k<MK;k++)
				printf("%d ", path2[i*MK*MK+j*MK+k]);
			printf("\n\t");
		}
		printf("\n");
	}
}
*/

	if(bst == 0)
	{	trace1.clear();
		trace1.resize(SSZ, vector<int>((int)breaks.size() / 2, -1));
		trace2 = trace1;
		if(ploidity == 1)
		{	double maxP;
			for(i = 0; i < MK; i++, lpp++)
			{	double mm = tP1[i];
				for(j = 1; j < MK; j++) mm = max(mm, tP1[j * MK + i]);
				*lpp = 0;
				for(j = 0; j < MK; j++) 
				{	cumsum[j] = exp(tP1[j * MK + i] - mm) * tpara.Vh[j];
					*lpp += cumsum[j];
				}
				*lpp = log(*lpp) + mm;
				if(i == 0 || maxP < *lpp) maxP = *lpp;

				k = _sample(cumsum, MK, tpara.maximum & (1<<1));
				int *ppp, es = i;
				if(stepn1 > 0) ppp = &path1[MK * MK * (stepn1 - 1)];
				for(j = bed - 2; j >= bst; j -= 2)
				{	trace1[i][j/2] = es;
					if(j > bst && breaks[j + 1] != 1) { es = *(ppp + k * MK + es); ppp -= MK * MK; }
				}
			}
			double *tlpp = lpp - SSZ;
			*lss = 0;
			for(i = 0; i < MK; i++, tlpp++) 
				*lss += exp((*tlpp) - maxP);
			*lss = log(*lss) + maxP;
		}
		else
		{	double P1[MK], P2[MK];
			int as1[MK], as2[MK];
			for(i = 0; i < MK; i++)
			{	double mm1 = tP1[i], mm2 = tP2[i];
				for(j = 1; j < MK; j++)
				{	mm1 = max(mm1, tP1[j * MK + i]);
					mm2 = max(mm2, tP2[j * MK + i]);
				}
				P1[i] = P2[i] = 0;
				double cumsum2[MK];
				for(j = 0; j < MK; j++)
				{	cumsum[j] = exp(tP1[j * MK + i] - mm1);
					cumsum2[j] = exp(tP2[j * MK + i] - mm2);
					P1[i] += cumsum[j];
					P2[i] += cumsum2[j];
				}
				P1[i] = log(P1[i]) + mm1 + log(tpara.Vh[i] + 1e-20);
				P2[i] = log(P2[i]) + mm2 + log(tpara.Vh[i] + 1e-20);

				as1[i] = _sample(cumsum, MK, tpara.maximum & (1<<1));
				as2[i] = _sample(cumsum2, MK, tpara.maximum & (1<<1));
			}
			double maxP;
			vector<vector<int> >::const_iterator iss = ss.begin();
			for(i = 0; i < SSZ; i++, lpp++, iss++)
			{	j = (*iss)[0]; k = (*iss)[1];
				double a = 0.5, b = 0.5;
				if(j == k) *lpp = P1[j] + P2[k];
				else
				{	double mm = max(P1[j] + P2[k], P1[k] + P2[j]);
					a = exp(P1[j]+P2[k]-mm); b = exp(P1[k]+P2[j]-mm);
					*lpp = log(a+b) + mm; 
				}
				if(i == 0 || maxP < *lpp) maxP = *lpp;

				int *ppp, es1 = j, es2 = k;
				if(gsl_ran_flat(gammar, 0, a + b) > a) { es1 = k; es2 = j; }
				int ts1 = as1[es1], ts2 = as2[es2];
				if(stepn1 > 0) ppp = &path1[MK * MK * (stepn1 - 1)];
				for(j = bed - 2; j >= bst; j -= 2)
				{	trace1[i][j/2] = es1;
					if(j > bst && breaks[j + 1] != 1) { es1 = *(ppp + ts1 * MK + es1); ppp -= MK * MK; }
				}
				if(stepn2 > 0) ppp = &path2[MK * MK * (stepn2 - 1)];
				for(j = bed - 2; j >= bst; j -= 2)
				{	trace2[i][j/2] = es2;
					if(j > bst && breaks[j + 1] != 0) { es2 = *(ppp + ts2 * MK + es2); ppp -= MK * MK; }
				}
			}
			double *tlpp = lpp - SSZ, *lmm1 = lmm, *lmm2, f;
			for(i = 0; i < MK; i++) *(lmm + i) = 0;
			*lss = 0;
			for(i = 0; i < MK; i++, lmm1++)
			{	lmm2 = lmm;
				for(j = 0; j <= i; j++, lmm2++, tlpp++)
				{	f = exp(*tlpp - maxP);
					(*lmm1) += f;
					(*lmm2) += f;
					(*lss) += f;
				}
			}
			for(i = 0; i < MK; i++) *(lmm + i) = log(*(lmm + i)) + maxP;
			(*lss) = log(*lss) + maxP;
		}
	}
	else
	{	double recp = *rp;
		double *tlpp = lpp - SSZ, *tlss = lss - 1;
		if(ploidity == 1)
		{	double maxP;
			for(i = 0; i < MK; i++, lpp++)
			{	tP1[i] += log(exp(*tlpp - *tlss) * (1. - recp) + recp * tpara.Vh[0]) + (*tlss);
				double mm = tP1[i];
				for(j = 1; j < MK; j++) 
				{	tP1[j * MK + i] += log(exp(*(tlpp + j) - *tlss) * (1. - recp) + recp * tpara.Vh[j]) + (*tlss);
					mm = max(mm, tP1[j * MK + i]);
				}
				*lpp = 0;
				for(j = 0; j < MK; j++) 
				{	cumsum[j] = exp(tP1[j * MK + i] - mm);
					*lpp += cumsum[j];
				}
				*lpp = log(*lpp) + mm;

				if(i == 0 || maxP < *lpp) maxP = *lpp;
			
				k = _sample(cumsum, MK, tpara.maximum & (1<<1));
				int *ppp, es = i;
				if(stepn1 > 0) ppp = &path1[MK * MK * (stepn1 - 1)];
				for(j = bed - 2; j >= bst; j -= 2)
				{	trace1[i][j/2] = es;
					if(j > bst && breaks[j + 1] != 1) { es = *(ppp + k * MK + es); ppp -= MK * MK; }
				}
			}
			tlpp = lpp - SSZ;
			*lss = 0;
			for(i = 0; i < MK; i++, tlpp++) 
				*lss += exp((*tlpp) - maxP);
			*lss = log(*lss) + maxP;
		}
		else
		{	
			double nr2 = (1. - recp) * (1. - recp), r2 = 2. * recp * recp, rnr = recp * (1. - recp);
			double *tlmm = lmm - MK, *otlpp = tlpp;
			double maxP;
			vector<vector<int> >::const_iterator iss = ss.begin();
			for(i = 0; i < SSZ; i++, lpp++, iss++)
			{	int s1 = (*iss)[0], s2 = (*iss)[1];
				
				tlpp = otlpp;
				vector<vector<int> >::const_iterator tss = ss.begin();
				double tprob[SSZ], tmax;
				for(j = 0; j < SSZ; j++, tlpp++, tss++)
				{	int ts1 = (*tss)[0], ts2 = (*tss)[1];
					double f1 = tP1[ts1*MK + s1] + tP2[ts2*MK + s2];
					double f2 = tP1[ts1*MK + s2] + tP2[ts2*MK + s1];
					double f3 = tP1[ts2*MK + s1] + tP2[ts1*MK + s2];
					double f4 = tP1[ts2*MK + s2] + tP2[ts1*MK + s1];
					double mm = max(f1, f2), nn = max(f3, f4);
					mm = max(mm, nn);
					double f = log((exp(f1 - mm) + exp(f2 - mm) + exp(f3 - mm) + exp(f4 - mm)) / (1. + (double)(s1==s2))) + mm;
					tprob[j] = log(exp(*tlpp - *tlss) * nr2 + rnr * (exp(*(tlmm + ts1) - *tlss) * tpara.Vh[ts2] + exp(*(tlmm + ts2) - *tlss) * tpara.Vh[ts1]) / (1. + (double)(ts1==ts2)) + r2 * Vh2[j]) + f + (*tlss);
					if(j == 0 || tmax < tprob[j]) tmax = tprob[j];
				}
				*lpp = 0;
				for(j = 0; j < SSZ; j++)
				{	cumsum[j] = exp(tprob[j] - tmax);
					*lpp += cumsum[j];
				}
				*lpp = log(*lpp) + tmax;
				if(i == 0 || maxP < *lpp) maxP = *lpp;

				//traceback
				k = _sample(cumsum, SSZ, tpara.maximum & (1<<1));
if(k >= SSZ) {for(j=0;j<SSZ;j++) printf("%f,",cumsum[j]);printf("tmax=%f\n",tmax);exit(0);}
				int ts1 = ss[k][0], ts2 = ss[k][1];
				cumsum[0] = tP1[ts1*MK + s1] + tP2[ts2*MK + s2];
				cumsum[1] = tP1[ts1*MK + s2] + tP2[ts2*MK + s1];
				cumsum[2] = tP1[ts2*MK + s1] + tP2[ts1*MK + s2];
				cumsum[3] = tP1[ts2*MK + s2] + tP2[ts1*MK + s1];	
				double mm = max(cumsum[0], cumsum[1]), nn = max(cumsum[2], cumsum[3]);
				mm = max(mm, nn);
				for(j = 0; j < 4; j++) cumsum[j] = exp(cumsum[j] - mm);
				k = _sample(cumsum, 4, tpara.maximum & (1<<1));
if(k>3) printf("!!!%f,%f,%f,%f\n", cumsum[0],cumsum[1],cumsum[2],cumsum[3]),exit(0);
				if(k > 1) { j = ts1; ts1 = ts2; ts2 = j; }
				int es1 = s1, es2 = s2;
				if(k == 1 || k == 3) { es1 = s2; es2 = s1; }
				
				int *ppp;
				if(stepn1 > 0) ppp = &path1[MK * MK * (stepn1 - 1)];
				for(j = bed - 2; j >= bst; j -= 2)
				{	trace1[i][j/2] = es1;
					if(j > bst && breaks[j + 1] != 1) { es1 = *(ppp + ts1 * MK + es1); ppp -= MK * MK; }
				}
				if(stepn2 > 0) ppp = &path2[MK * MK * (stepn2 - 1)];
				for(j = bed - 2; j >= bst; j -= 2)
				{	trace2[i][j/2] = es2;
					if(j > bst && breaks[j + 1] != 0) { es2 = *(ppp + ts2 * MK + es2); ppp -= MK * MK; }
				}
			}
			tlpp = lpp - SSZ; 
			double *lmm1 = lmm, *lmm2, f;
			for(i = 0; i < MK; i++) *(lmm + i) = 0;
			*lss = 0;
			for(i = 0; i < MK; i++, lmm1++)
			{	lmm2 = lmm;
				for(j = 0; j <= i; j++, lmm2++, tlpp++)
				{	f = exp(*tlpp - maxP);
					(*lmm1) += f;
					(*lmm2) += f;
					(*lss) += f;
				}
			}
			for(i = 0; i < MK; i++) *(lmm + i) = log(*(lmm + i)) + maxP;
			(*lss) = log(*lss) + maxP;
		}
	}
/*
if(id==12 && bst==0)
{	printf("trace1\n");
	for(i=0;i<SSZ;i++) 
	{	printf("%d,%d: ", ss[i][0], ss[i][1]);
		for(j = 0; j < (bed-bst)/2; j++) printf("%d ", trace1[i][j]);
		printf("\n");
	}
	printf("\ntrace2\n");
	for(i=0;i<SSZ;i++) 
	{	printf("%d,%d: ", ss[i][0], ss[i][1]);
		for(j = 0; j < (bed-bst)/2; j++) printf("%d ", trace2[i][j]);
		printf("\n");
	}
	exit(0);
}*/

if(false)
{	printf("[%d:%d-%d] ", id, breaks[bst], breaks[bed]);
	for(i = 0; i < SSZ; i++) printf("(%d,%d):%f ", ss[i][0], ss[i][1], *(lpp-SSZ+i));
	printf(" | lss=%f\n", *lss);
}
//	printf("lss=%f\n", *lss);fflush(stdout);

	lmm += MK;
	lss ++;
	k = breaks[bed] - breaks[bst];
	ipm += k;
	nhp += k * tpara.maxK;
	aszp += k;
	rp += k;
	dp += k;
	if(ploidity == 2) dp2 += k;
}

//strand = 0 or 1
int popTensor::_oneStrandProb(int id, vector<int> const &breaks, int bst, int bed, int strand, TENSORPARA const &tpara, double const *rp, int const *nhp, vector<vector<vector<float> > >::const_iterator ipm, int const *aszp, int MK, double *tP, double *sP, int *path)
{	int i, j, k, l, m;
	double *ddd = dataP + id * dSZ;
	int dist[(bed - bst) / 2];
	for(i = bst; i < bed; i += 2)
	{	dist[(i - bst) / 2] = breaks[i] - breaks[bst];
	}
	int stepn = 0;
	int *ppp = path;
	double tmpP[MK * MK];
	for(i = bst; i < bed; i += 2)
	{	for(j = i + 2; j < bed; j+=2)
		{	if(breaks[j+1]!=1-strand) break;
		}
		for(k = 0; k < MK; k++)
		{	double maxlp = -100000000.;
			double prob[(*aszp)];
			for(l = 0; l < (*aszp); l++)
			{	prob[l] = 0;
				for(m = i; m < j; m += 2)
				{	int d = dist[(m-bst)/2];
					prob[l] += _getLp1(ddd + m / 2 * (*aszp), k, l, (m==i), tpara.logscale, ipm + d, nhp + d * tpara.maxK, aszp + d, tpara.maxK, tpara.probAlpha);
				}
				maxlp = max(maxlp, prob[l]);
			}
			double f = 0;
			for(l = 0; l < (*aszp); l++) f += exp(prob[l] - maxlp);
			f = log(f) + maxlp;

			if(i == bst) 
			{	for(l = 0; l < MK; l++) tP[k * MK + l] = -100000000.;
				tP[k * MK + k] = /*log(tpara.Vh[k]) + */f;
			}
			else
			{	double recp = *(rp + dist[(i-bst)/2]), recpv = recp * tpara.Vh[k];
				for(l = 0; l < MK; l++)
				{	if(recpv>0) tP[l * MK + k] = log(recpv + (1. - recp) * exp(tP[l * MK + k] - sP[l])) + sP[l] + f;
					else tP[l * MK + k] += log(1. - recp) + f;
if(tP[l*MK+k]>-1000000000 && tP[l*MK+k]<1000000) ;
else 
{printf("%d:%d-%d!!!",id,i,j);
printf("MK=%d, recp=%f, f=%f, sP[%d]=%f\n", MK, recp, f, l, sP[l]);
for(l=0;l<(*aszp);l++) printf("%f,",prob[l]);
printf("|%f\n",maxlp);
exit(0);
}
					double cumsum[MK];
					for(m = 0; m < MK; m++) 
					{	if(m != k) cumsum[m] = tmpP[l * MK + m];
						else cumsum[m] = tmpP[l * MK + m] *  (1. + (1. - recp) / recpv);
					}
					*(ppp + l * MK + k) = _sample(cumsum, MK, tpara.maximum & (1<<1));
				}
				if(k==MK-1) 
				{	stepn++;
					ppp += MK * MK;
				}
			}
		}
		double *tp1 = &tP[0], *tp2 = &tmpP[0];
		for(k = 0; k < MK; k++)
		{	double m = *tp1;
			for(l = 1; l < MK; l++) m = max(m, *(tp1 + l));
			sP[k] = 0;
			for(l = 0; l < MK; l++, tp1++, tp2++) 
			{	*tp2 = exp(*tp1 - m);
				sP[k] += *tp2;
			}
			sP[k] = log(sP[k]) + m;
		}

		i = j - 2;
	}

	return stepn;
}

int popTensor::_tracebackOne(int id, vector<int> const &breaks, int bst, int bed, TENSORPARA &tpara, int ploidity, int *&pp, int *&pp2, bool *&rr, bool *&rr2, double const *&rp, int const *&nhp, vector<vector<vector<float> > >::const_iterator &ipm, double *&lpp, double *&lmm, double *&lss, int SSZ, int MK, vector<vector<int> > const &ss, vector<vector<int> > const &mlist, double *Vh2)
{
	//trace back needs to update rec
	int i, j, k, flip = 0;
	int s1 = min(*pp, tpara.maxK), s2, code = s1; 
	if(ploidity == 2)
	{	s2 = min(*pp2, tpara.maxK);
		code = s1 * (s1 + 1) / 2 + s2;
		if(s1 < s2) code = s2 * (s2 + 1) / 2 + s1;
		if(s1 == s2) flip = (int)gsl_ran_flat(gammar, 0., 2.);
		else if(trace1[code][bed/2-1] == s2 && trace2[code][bed/2-1] == s1) flip = 1;
	}

	double const *orp = rp;
	vector<int>::iterator itr;
	if(flip == 0) itr = trace1[code].begin();
	else itr = trace2[code].begin();
	for(i = bed - 2; i >= bst; i -= 2)
	{	int st = breaks[i], ed = breaks[i + 2], stat = *(itr + i / 2), ts;

                if(*rr || i == (int)breaks.size() - 4)
		{	ts = stat;
			if(ts >= tpara.maxK && (tpara.maximum & (1<<1))==0)
			{	do { ts = tpara.maxK + (int)(- log(gsl_ran_flat(gammar, 0, tpara.Vh[tpara.maxK]) / tpara.Vh[tpara.maxK]) / log((1. + tpara.A) / tpara.A));
				} while(tpara.maxHapK > 0 && ts >= tpara.maxHapK && ts > tpara.maxK);
			}
		}
		else ts = *(pp+1);
		rr--;
		for(j = st; j < ed; j++, pp--,rr--,rp--)
		{	*pp = ts; *rr = false; }
		rr++;
		if(i > bst)
		{	if(*(itr + i / 2 - 1) != stat) *rr = true;
			else if(breaks[i + 1] != 1-flip)
			{	double f = (*rp) * tpara.Vh[stat] / ((*rp) * tpara.Vh[stat] + 1. - (*rp));
				if((tpara.maximum & (1<<2))!=0 && f > 0.5) *rr = true;
				else if((tpara.maximum & (1<<2))==0 && gsl_ran_flat(gammar, 0, 1.) <= f) *rr = true;	
			}
		}
	}

	if(ploidity == 2)
	{	rp = orp;
		itr = trace2[code].begin();
		if(flip == 0) itr = trace2[code].begin();
		else itr = trace1[code].begin();
		for(i = bed - 2; i >= bst; i -= 2)
		{	int st = breaks[i], ed = breaks[i + 2], stat = *(itr + i / 2), ts;
	                if(*rr2 || i == (int)breaks.size() - 4)
			{	ts = stat;
				if(ts >= tpara.maxK && (tpara.maximum & (1<<1))==0)
				{	do { ts = tpara.maxK + (int)(- log(gsl_ran_flat(gammar, 0, tpara.Vh[tpara.maxK]) / tpara.Vh[tpara.maxK]) / log((1. + tpara.A) / tpara.A));
					} while(tpara.maxHapK > 0 && ts >= tpara.maxHapK && ts > tpara.maxK);
				}
			}
			else ts = *(pp2+1);
			rr2--;
			for(j = st; j < ed; j++, pp2--,rr2--,rp--)
			{	*pp2 = ts; *rr2 = false; }
			rr2++;
			if(i > bst)
			{	if(*(itr + i / 2 - 1) != stat) *rr2 = true;
				else if(breaks[i + 1] != flip)
				{	double f = (*rp) * tpara.Vh[stat] / ((*rp) * tpara.Vh[stat] + 1. - (*rp));
					if((tpara.maximum & (1<<2))!=0 && f > 0.5) *rr2 = true;
					else if((tpara.maximum & (1<<2))==0 && gsl_ran_flat(gammar, 0, 1.) <= f) *rr2 = true;	
				}
			}
		}
	}

	k = breaks[bed] - breaks[bst];
	nhp -= k*tpara.maxK;
	ipm -= k;
	lpp -= SSZ;
	lmm -= MK;
	lss --;

	return flip;
}

double popTensor::_getLp(int id, vector<int> const &breaks, int bid, bool flip, int ploidity, TENSORPARA const &tpara, float const *dp, float const *dp2, int const *pp, int const *pp2, vector<vector<vector<float> > >::const_iterator ipm, int const *nhp, int const *aszp)
{
	double a = tpara.probAlpha, prob;
	
	int p1 = min(tpara.maxK, (*pp)), p2; 
	if(ploidity == 2) p2 = min(tpara.maxK, (*pp2));
	
	double *dddp = dataP + id * dSZ + bid / 2 * (*aszp);
	if(flip == 1) dddp += dSZ;
	prob = _getLp1(dddp, p1, (int)(*dp), (breaks[bid+1]!=1-flip), tpara.logscale, ipm, nhp, aszp, tpara.maxK, a);

	if(ploidity == 2)
	{	if(flip == 1) dddp -= dSZ;
		else dddp += dSZ;
		double prob2 = _getLp1(dddp, p2, (int)(*dp2), (breaks[bid+1]!=flip), tpara.logscale, ipm, nhp, aszp, tpara.maxK, a);
		if(tpara.logscale) prob += prob2;
		else prob *= prob2;
	}	
	
	return prob;
}

double popTensor::_getLp1(double *dddp, int pop, int state, bool add, bool logscale, vector<vector<vector<float> > >::const_iterator ipm, int const *nhp, int const *aszp, int maxK, double a)
{
	double f, t;
	f = *(dddp + state);
	if(add)
	{	if(pop >= maxK) t = 1. / (double)(*aszp);
		else t = ((double)(*ipm)[pop][state] + a) / ((double)(*(nhp + pop)) + a * (double)(*aszp));
		if(logscale) f += log(t);
		else f *= t;
	}

	return f;
}

double popTensor::_imputeOne(int id, vector<int> const &breaks, int bst, int bed, int flip, int ploidity, TENSORPARA const &tpara, float *dp, float *dp2, int const *pp, int const *pp2, vector<vector<vector<float> > >::const_iterator ipm, int const *nhp, int const *aszp)
{
	int i, j, k, l;
	int dist[(bed - bst) / 2];
	for(i = bst; i < bed; i += 2)
	{	dist[(i - bst) / 2] = breaks[i] - breaks[bst];
	}
	double *ddd = dataP + id * dSZ;
	if(flip == 1) ddd += dSZ;

	double RT = 0;
	for(i = bst; i < bed; i += 2)
	{	for(j = i + 2; j < bed; j+=2)
		{	if(breaks[j+1]!=1-flip) break;
		}
		double prob[(*aszp)], mpp = -100000000.;
		for(l = 0; l < (*aszp); l++)
		{	prob[l] = 0;
			for(k = i; k < j; k += 2)
			{	int d = dist[(k-bst)/2];
				prob[l] += _getLp1(ddd + k / 2 * (*aszp), *pp, l, (k==i), tpara.logscale, ipm + d, nhp + d * tpara.maxK, aszp + d, tpara.maxK, tpara.probAlpha);
			}
			mpp = max(mpp, prob[l]);
		}
		for(l = 0; l < (*aszp); l++) prob[l] = exp(prob[l] - mpp);
		k = _sample(prob, (*aszp), tpara.maximum & 1);
		RT += mpp;
if(k >= (*aszp))
{	for(l = 0; l < (*aszp); l++) printf("%d:%f ", l, prob[l]);
	printf("| k=%d, i=%d\n", k, i);
	exit(0);
}
		int st = breaks[i], ed = breaks[j];
		for(l = st; l < ed; l++, dp++) *dp = k;
		i = j - 2;
	}
	if(ploidity == 2)
	{	if(flip == 1) ddd -= dSZ;
		else ddd += dSZ;
		for(i = bst; i < bed; i += 2)
		{	for(j = i + 2; j < bed; j+=2)
			{	if(breaks[j+1]!=flip) break;
			}
			double prob[(*aszp)], mpp = -100000000.;
			for(l = 0; l < (*aszp); l++)
			{	prob[l] = 0;
				for(k = i; k < j; k += 2)
				{	int d = dist[(k-bst)/2];
					prob[l] += _getLp1(ddd + k / 2 * (*aszp), *pp2, l, (k==i), tpara.logscale, ipm + d, nhp + d * tpara.maxK, aszp + d, tpara.maxK, tpara.probAlpha);
				}
				mpp = max(mpp, prob[l]);
			}
			for(l = 0; l < (*aszp); l++) prob[l] = exp(prob[l] - mpp);
			k = _sample(prob, (*aszp), tpara.maximum & 1);
			RT += mpp;
			int st = breaks[i], ed = breaks[j];
			for(l = st; l < ed; l++, dp2++) *dp2 = k;
			i = j - 2;
		}
	}
	return(RT);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void popTensor::readInput(MYDATA &mydata, char const *input, char const *fS, vector<double> &recomb)
{	int i, j, k, l;
	char fname[200], tmp[1000];
	sprintf(fname, "%s.popr", input);
	FILE *f = fopen(fname, "r");
	recomb.clear();
	i = 0;
	while(fgets(tmp, 1000, f) != NULL)
	{	if((int)strlen(tmp) < 3) continue;
		for(j = 0; j < (int)strlen(tmp); j++) if(tmp[j] == ' ' || tmp[j] == '\t') break;
		recomb.push_back(atof(&tmp[j + 1]));
	}
	fclose(f);
	mydata.L = (int)recomb.size();

	string str;
	ifstream inFile;
	sprintf(fname, "%s.pop", input);
        inFile.open(fname);
	vector<vector<int> > stateShort;
	vector<int> row;
	int id, last = -1, counter = 0;
	while(true)
	{	inFile >> str;
		if(inFile.eof()) break;
		if(str[0] == 'i' && str[1] == 'd')
		{	if(last >= 0)
			{	if(stateShort.size() <= id) stateShort.push_back(row);
				else stateShort[id] = row;
			}
			id = atoi(&str[2]);
			if(id == 0 && last >= 0) 
			{	counter++;
			//	if(counter == 40) break;
			}
			last = 0;
			row.clear();
		}
		else
		{	for(i = 0; i < str.size(); i++) if(str[i] == ':') break;
			int pos = atoi(&str[0]);
			row.push_back(pos);
			row.push_back(atoi(&str[i + 1]));
			last = pos;
		}
	}
	inFile.close();

	mydata.data = NULL;
	if(last >= 0)
	{	if((int)stateShort.size() <= id) stateShort.push_back(row);
                else stateShort[id] = row;

                vector<double> logp;
		mydata.indN = (int)stateShort.size();
                
		mydata.data = new float[mydata.indN * mydata.L];
		int maxK = 0;
		float *p = mydata.data;
                for(i = 0; i < mydata.indN; i++) 
                {       for(j = 0; j < (int)stateShort[i].size(); j+=2) 
			{	int st = stateShort[i][j], ed = mydata.L, s = stateShort[i][j + 1];
				if(j < (int)stateShort[i].size() - 2) ed = stateShort[i][j + 2];
				for(k = st; k < ed; k++, p++) (*p) = s;
				maxK = max(maxK, s + 1);
			}
                }

		dSZ = 0;
		mydata.breaks.clear();
        	mydata.breaks.resize(mydata.indN / mydata.ploidity);
                for(i = 0; i < (int)mydata.breaks.size(); i++)
                {       int j1 = 0, j2 = 0;
			vector<int>::iterator ist1 = stateShort[i*mydata.ploidity].begin(), ist2 = stateShort[i*mydata.ploidity+mydata.ploidity - 1].begin();
			int ssz1 = stateShort[i*mydata.ploidity].size(), ssz2 = stateShort[i*mydata.ploidity+mydata.ploidity-1].size();
			while(j1 < ssz1 && j2 < ssz2)
			{	int to = - 1;
				if(j1 >= ssz1) to = 1;
				else if(j2 >= ssz2) to = 0;
				else if(*(ist1+j1) < *(ist2+j2)) to = 0;
				else if(*(ist1+j1) > *(ist2+j2)) to = 1;
				else to = 2;
				if(to == 0 || to == 2) { mydata.breaks[i].push_back(*(ist1+j1)); j1 += 2; }
				else if(to == 1) { mydata.breaks[i].push_back(*(ist2+j2)); j2 += 2; }
				if(to == 2) j2 += 2;
				mydata.breaks[i].push_back(to);
			}
                        mydata.breaks[i].push_back(mydata.L);
                        mydata.breaks[i].push_back(2);
			dSZ = max(dSZ, (int)mydata.breaks[i].size() / 2 - 1);
		}
		dSZ *= maxK;

		char *haps = new char[mydata.indN * mydata.L];
		_readHaplotype(input, haps, mydata.indN, mydata.L);	
		vector<vector<int> > ann(mydata.L, vector<int>(maxK, 0)), nnn = ann;
		for(i = 0; i < mydata.L; i++)
		{	vector<int>::iterator ian = ann[i].begin(), inn = nnn[i].begin();
			char *happ = haps + i * mydata.indN;
			for(j = 0; j < mydata.indN; j++, happ++)
			{	int s = (int)mydata.data[j * mydata.L + i];
				*(ian + s) = *(ian + s) + (int)(*happ);
				*(inn + s) = *(inn + s) + 1;
			}
		}
		if(dataP != NULL) delete dataP;
		dataP = new double[mydata.indN * dSZ];
		for(i = 0; i < mydata.indN; i++)
		{	int lastpos = 0, bsz = (int)mydata.breaks[i / mydata.ploidity].size();
			double *lp = dataP + dSZ * i;
			for(j = 0; j < bsz - 2; j+=2)
			{	int st = mydata.breaks[i / mydata.ploidity][j];
				int ed = mydata.breaks[i / mydata.ploidity][j + 2];
				for(k = 0; k < maxK; k++)
				{	double s = 0;
					for(l = lastpos; l < ed; l++)
                                	{       double fq = ((double)ann[l][k] + 0.5) / ((double)nnn[l][k] + 1.);
						double hh = (double)haps[l*mydata.indN+i];
                                        	if(mydata.ploidity == 2) s += log(hh * fq + (1. - hh) * (1. - fq));
                                		else if((int)hh == 0) s += log(1. - fq) * 2.;
						else if((int)hh == 1) s += log(fq) + log(1. - fq);
						else s += log(fq) * 2.;
if(s > -10000 || s < 10000) ;
else
{	printf("i=%d,j=%d,k=%d, s=%f, ploidity=%d, hh=%d, fq=%f, haps=%d, ann=%d, nnn=%d\n", i, j, k, s, mydata.ploidity, (int)hh, fq, haps[l*mydata.indN+i],ann[l][k],nnn[l][k]);
}
					}
					lp[(j/2) * maxK + k] = s;
if(s > -10000 || s < 10000) ;
else
{	printf("i=%d,j=%d,k=%d, s=%f\n", i, j, k, s);
exit(0);
}
				}
				lastpos = ed;
			}
		}
		delete haps;

		mydata.fixAllele.clear();
		mydata.fixState.clear();
		if(fS != NULL) 
		{	_readPop(fS, mydata.fixState);
			if((int)mydata.fixState.size() < mydata.indN) mydata.fixState.resize(mydata.indN, vector<int>());
		}
	}
	mydata.indIndex = new int[mydata.indN/mydata.ploidity + 1];
	for(i = 0; i < mydata.indN/mydata.ploidity + 1; i++) mydata.indIndex[i] = i;
	mydata.totalN = mydata.indN/mydata.ploidity;
}

void popTensor::_readHaplotype(char const *fname, char *haps, int N, int L)
{	int i, j, k;
	char *p = haps;
	char str[1000];
	FILE *f;
	vector<char> aa(2, 0);
	aa[1] = 1;
	vector<vector<char> > alleles(L, aa);
	sprintf(str, "%s.snps", fname);
	f = fopen(str, "r");
	for(i = 0; i < L; i++)
	{	char tmp[1000];
		fgets(tmp, 1000, f);
		k=0;
		int l = (int)strlen(tmp) - 2;
		for(j = 0; j < l; j++)
		{	if(tmp[j] == ' ' || tmp[j] == '\t')
			{	k++;
				if(k == 4) break;
			}
		}
		if(j < l) 
		{	alleles[i][0] = tmp[j + 1];
			alleles[i][1] = tmp[j + 3];
		}
	}	
	fclose(f);

	sprintf(str, "%s.g", fname);
	f = fopen(str, "r");
	char c, tmp1[100], tmp2[100], tmp3[100];
	fscanf(f, "%s %s %s ", tmp1, tmp2, tmp3);
	for(i = 0; i < N; i++) fscanf(f, "%s ", tmp1, tmp2);
	fscanf(f, "\n");
	for(i = 0; i < L; i++)
	{	fscanf(f, "%s %s %d ", tmp1, tmp2, &k);
		for(j = 0; j < N; j++)
		{	fscanf(f, "%c ", &c);
			if(c == alleles[i][0]) *p = 0;
			else if(c == alleles[i][1]) *p = 1;
			else printf("!!! alleles do not match at position %d: %c != %c or %c\n", i, c, alleles[i][0], alleles[i][1]),exit(0);
			p++;
		}
		fscanf(f, "\n");
	}
	fclose(f);
}
