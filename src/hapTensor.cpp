#include "hapTensor.h"
#define STDFORMAT 1 

hapTensor::hapTensor(unsigned int rs):tensorHMMbase(rs)
{
	dataP = nullP = errP = NULL;
}

hapTensor::~hapTensor()
{	if(dataP != NULL) delete dataP;
	if(nullP != NULL) delete nullP;
	if(errP != NULL) delete errP;
}

void hapTensor::_getAsz(MYDATA const &mydata)
{
	int i;
        for(i = 0; i < mydata.L; i++) mydata.asz[i] = 2;
}

/////////////////////////////////////////////////////////////////////////////////////
void hapTensor::run(char const *finput, bool likelihood, int burnin, int mcmc, int maxHapK, int ploidity, double err, double AA, bool samplemax, int cut, char const *foutput, vector<string> const &fS, vector<string> const &fA, bool sqc)
{
	MYDATA mydata;
	
	trueploidity = mydata.ploidity = ploidity;
	mydata.L = -1;
	double lambda = readInput(mydata, finput, fS, fA, err, cut, likelihood); //need revision for fix state and fix alleles
	
	TENSORPARA tpara;
	initPara(tpara);
	tpara.probAlpha = min(1., max(1e-1, 10./ (double)mydata.indN / lambda));
	tpara.maxHapK = maxHapK;
	tpara.burnin = burnin;
	tpara.mcmc = mcmc;
	tpara.samplemaximum = samplemax;
	tpara.samplefull = true;
	tpara.changeAlpha = 1;
	tpara.logscale = false;
	if(AA > 0) tpara.probAlpha = AA;

	inferStructure(mydata, tpara, foutput, sqc);
	
	char fname[200];
	if(tpara.singlecall) sprintf(fname, "%s.g0", foutput);
	else sprintf(fname, "%s.g", foutput);
	outputHaplotype(fname, mydata);
	
	if(mydata.data != NULL) delete mydata.data;
	if(mydata.indIndex != NULL) delete mydata.indIndex;
}

/////////////////////////////////////////////////////////////////////////////////////
void hapTensor::_updatePara_remove(vector<float> &para, float data, float b, int id, int pos)
{	para[(int)data] -= b;	
}

void hapTensor::_updatePara_add(vector<float> &para, float data, float b, int id, int pos)
{	para[(int)data] += b;
}

//need to update pointers, bst and bed are breaks index, if they exist, otherwise they are positions
void hapTensor::_sumforwardOne(int id, vector<int> const &breaks, int bst, int bed, TENSORPARA &tpara, int ploidity, float *&dp, float *&dp2, double const *&rp, int const *&nhp, vector<vector<vector<float> > >::const_iterator &ipm, int const *&aszp, double *&lpp, double *&lmm, double *&lss, int SSZ, int MK, vector<vector<int> > const &ss, vector<vector<int> > const &mlist, double *Vh2)
{	int i, j;
	double f, t;	

	(*lss) = 0;
	if(ploidity == 2)
	{	for(j = 0; j < MK; j++) *(lmm + j) = 0;	}

	if(bst <= 0)
	{	for(j = 0; j < SSZ; j++, lpp++)
        	{       (*lpp) = 0;
			f = _getLp(id, breaks, bst, false, ploidity, tpara, dp, dp2, &ss[j][0], &ss[j][1], ipm, nhp, aszp);
			if(ploidity == 1) t = tpara.Vh[ss[j][0]];
			else t = Vh2[j] * (1. + (double)(ss[j][0] != ss[j][1]));
			(*lpp) = t * f;
			(*lss) += (*lpp);
			if(ploidity == 2)
			{	*(lmm + ss[j][0]) += *lpp;
				*(lmm + ss[j][1]) += *lpp;
			}
        	}
	}
	else
	{	double recp = *rp;
		double recp2 = recp * recp;
		double nrecp2 = (1. - recp) * (1. - recp);
		double *tmp;
		if(tpara.singlecall) { recp = recp2 = 1.-1e-10;nrecp2 = 1e-10; }
		
                double *tlpp = lpp - SSZ;
		if(ploidity == 1)
		{	for(j = 0; j < SSZ; j++, lpp++, tlpp++)
			{	f = _getLp(id, breaks, bst, false, ploidity, tpara, dp, dp2, &ss[j][0], &ss[j][1], ipm, nhp, aszp);
        	                t = tpara.Vh[j] * (*(lss - 1)) * recp + (*tlpp) * (1. - recp);
	                        (*lpp) = t * f;
        	                (*lss) += (*lpp);
	                }
		}
		else
		{	double rrs = recp2 * (*(lss - 1)), rnr = recp * (1. - recp), nr2 = nrecp2; 
			int s1, s2;
			rrs /= nr2; rnr /= nr2;
	
			double tppp[3] = {1., 1., 1.};
			double *ppp = dataP + (id / trueploidity * dSZ + bst) * 3, *nnp = nullP + (id / trueploidity * dSZ + bst) * 2, ep = errP[bst];
			double fq[MK], pf[MK], tf[MK];
			double aa = tpara.defmut * tpara.probAlpha, nn = tpara.probAlpha;
			tmp = lmm - MK;
			for(j = 0; j < MK; j++, tmp++)
			{	fq[j] = tpara.defmut;
				if(j < tpara.maxK && (double)(*(nhp + j)) + nn > 0) fq[j] = ((double)(*ipm)[j][1] + aa) / ((double)*(nhp + j) + nn);
				pf[j] = (*ppp) / 2. + fq[j] * (- (*ppp) + *(ppp+1));
				pf[j] *= nr2;
				tf[j] = (*tmp) * rnr;
				if(tpara.Vh[j] > 0) tf[j] /= (tpara.Vh[j]);
			}
			double psum = (*ppp) - 2. * (*(ppp+1)) + (*(ppp+2));
			double epnp = ep * (*nnp + (*(nnp + 1))), p[3], t;
			epnp*=nr2;
			for(j = 0; j < MK; j++) pf[j] *= (1. - ep);
			psum *= (1. - ep);
			psum*=nr2;

			vector<vector<int> >::const_iterator iss = ss.begin();
			double *ivh2 = Vh2;
			double *ppf1 = &pf[0], *qqf1 = &fq[0], *ttf1 = &tf[0], *mmp1 = lmm, rrs2 = 2. * rrs;
			for(s1 = 0; s1 < MK; s1++, ppf1++, qqf1++, ttf1++, mmp1++)
			{	double *ppf2 = &pf[0], *qqf2 = &fq[0], *ttf2 = &tf[0], *mmp2 = lmm, t0 = (*ppf1) + epnp, tt = (*qqf1) * psum, ttt = rrs2 + (*ttf1);
				for(s2 = 0; s2 < s1; s2++, ppf2++, qqf2++, ttf2++, mmp2++, lpp++, tlpp++, ivh2++)
				{	f = t0 + *ppf2 + tt * (*qqf2);
					t = (ttt + (*ttf2)) * (*ivh2);
					t += (*tlpp);
					(*lpp) = t*f;
					(*mmp1) += (*lpp);
					(*mmp2) += (*lpp);
				}
				f = t0 + *ppf2 + tt * (*qqf2);
				t = (rrs + (*ttf1)) * (*ivh2);
                                t += (*tlpp);
                                (*lpp) = t*f;
                                (*mmp1) += (*lpp);
                                (*mmp2) += (*lpp);
				lpp++; tlpp++; ivh2++;
			}
			mmp1=lmm;
			for(j = 0; j < MK; j++, mmp1++) (*lss) += (*mmp1);
			(*lss) /= 2.;
		}
                if((*lss) < 1e-10)
		{       tlpp = lpp - SSZ;
                        for(j = 0; j < SSZ; j++, tlpp++) (*tlpp) *= 1e+20;
			if(ploidity == 2)
			{	tmp = lmm;
				for(j = 0; j < MK; j++, tmp++)
					(*tmp) *= 1e+20;
			}
			(*lss) *= 1e+20;
		}
	}

	lmm += MK;
	lss++;
	ipm++;
	nhp += tpara.maxK;
	aszp++;
	rp++;
	dp++;
	if(ploidity == 2) dp2++;
}

int hapTensor::_tracebackOne(int id, vector<int> const &breaks, int bst, int bed, TENSORPARA &tpara, int ploidity, int *&pp, int *&pp2, bool *&rr, bool *&rr2, double const *&rp, int const *&nhp, vector<vector<vector<float> > >::const_iterator &ipm, double *&lpp, double *&lmm, double *&lss, int SSZ, int MK, vector<vector<int> > const &ss, vector<vector<int> > const &mlist, double *Vh2)
{
	pp--;
	nhp -= tpara.maxK;
	ipm--;
	lpp -= SSZ;
	lmm -= MK;
	lss --;
	rr--;
	rp--; 
	if(ploidity == 2) 
	{	pp2 --;
		rr2--;
	}
	return 0;
}

double hapTensor::_getLp(int id, vector<int> const &breaks, int bid, bool flip, int ploidity, TENSORPARA const &tpara, float const *dp, float const *dp2, int const *pp, int const *pp2, vector<vector<vector<float> > >::const_iterator ipm, int const *nhp, int const *aszp)
{
	double f1, f2, p[3];
	double a = tpara.defmut * tpara.probAlpha, n = tpara.probAlpha;
        f1 = f2 = a / n;
		
        if(*pp < tpara.maxK)
        {       f1 = ((double)(*ipm)[*pp][1] + a) / ((double)*(nhp + (*pp)) * (3. - (double)ploidity) + n);
                if((double)(*(nhp + (*pp))) + n == 0) f1 = tpara.defmut;
        }
        f2 = f1;
        if(ploidity == 2 && (*pp2) < tpara.maxK)
        {       f2 = ((double)(*ipm)[*pp2][1] + a) / ((double)*(nhp + (*pp2)) + n);
                if((double)(*(nhp + (*pp2))) + n == 0) f2 = tpara.defmut;
        }
	if(trueploidity > ploidity)
	{	p[0] = 1. - f1;
		p[1] = f1;
		p[2] = 0;
	}
	else
	{	p[0] = (1. - f1) * (1. - f2);
	        p[1] = f1 * (1. - f2) + (1. - f1) * f2;
	        p[2] = f1 * f2;
	}

	double prob;
	if(dp != NULL && (!tpara.imputedata || trueploidity > ploidity)) prob = p[(int)(*dp)];	
	else
	{	double *ppp = dataP + ((int)(id/trueploidity)*dSZ + bid) * 3, *nnp = nullP + ((int)(id/trueploidity)*dSZ + bid) * 2, ep = errP[bid];
		prob = ep * (*nnp + (*(nnp+1)));
		prob += (1. - ep) * (p[0] * (*ppp) + p[1] * (*(ppp + 1)) + p[2] * (*(ppp + 2)));
	}
	if(prob<=0 || prob > 1.+1e-10) printf("!!!id=%d, pos=%d, prob=%f\n", id, bid, prob),exit(0);
	if(tpara.logscale) return(log(prob));
	else return(prob);
}

double hapTensor::_imputeOne(int id, vector<int> const &breaks, int bst, int bed, int flip, int ploidity, TENSORPARA const &tpara, float *dp, float *dp2, int const *pp, int const *pp2, vector<vector<vector<float> > >::const_iterator ipm, int const *nhp, int const *aszp)
{
	int i;
	double f1, f2, p[3];
	double a = tpara.defmut * tpara.probAlpha, n = tpara.probAlpha;
        f1 = f2 = a / n;
		
        if(*pp < tpara.maxK)
        {       f1 = ((double)(*ipm)[*pp][1] + a) / ((double)*(nhp + (*pp)) * (3. - (double)ploidity) + n);
                if((double)(*(nhp + (*pp))) + n == 0) f1 = tpara.defmut;
        }
        f2 = f1;
        if(ploidity == 2 && (*pp2) < tpara.maxK)
        {       f2 = ((double)(*ipm)[*pp2][1] + a) / ((double)*(nhp + (*pp2)) + n);
                if((double)(*(nhp + (*pp2))) + n == 0) f2 = tpara.defmut;
        }
	p[0] = (1. - f1) * (1. - f2);
	p[1] = f1 * (1. - f2) + (1. - f1) * f2;
        p[2] = f1 * f2;
	
	double gsum[3], *ppp = dataP + (id / trueploidity * dSZ + bst) * 3, *nnp = nullP + (id / trueploidity * dSZ + bst) * 2, ep = errP[bst];
	for(i = 0; i < 3; i++)
	{	gsum[i] = (1. - ep) * p[i] * (*(ppp + i));
		if(i == 0) gsum[i] += ep * (*nnp);
		else if(i == 2) gsum[i] += ep * (*(nnp + 1));
	}
	i = _sample(gsum, 3, tpara.maximum & 1);

        if(ploidity == 1) (*dp) = i;
	else
	{	if(i == 0 || i == 2) (*dp) = (*dp2) = (int)(i / 2);
		else
		{	(*dp) = (*dp2) = 0;
			double un = gsl_ran_flat(gammar, 0, p[1]);
			if(*pp == *pp2)
			{	if(gsl_ran_flat(gammar, 0, 1.) < 0.5) (*dp) = 1;
				else (*dp2) = 1;
			}	
			else if(un <= f1 * (1. - f2) || ((tpara.maximum & 1)  && f1 * (1. - f2) > f2 * (1. - f1))) (*dp) = 1;
			else (*dp2) = 1;
		}
	}
	return(gsum[2]);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double hapTensor::readInput(MYDATA &mydata, char const *input, vector<string> const &fS, vector<string> const &fA, double error, int cut, bool likelihood)
{
	char fname[200];
	sprintf(fname, "%s.data", input); 
	if(STDFORMAT == 1) sprintf(fname, "%s.snps", input);
	FILE *f;
	char tmp[10000];
	int i, j, k;
	
	//read snp info
	f = fopen(fname, "r");
	if(f == NULL) printf("!!!cannot open file %s\n", fname), exit(0);
	mydata.snpinfo.clear();
	i = 0;
	while(fgets(tmp, 10000, f) != NULL)
	{	if(mydata.L < 0 || i < mydata.L)
		{	SNPINFO tsnp;
			for(j = 0; j < (int)strlen(tmp); j++) if(tmp[j] == ' ' || tmp[j] == '\t') break;
			tmp[j] = 0;
			tsnp.snpid = tmp;
			tmp[j] = ' ';
			tsnp.chr = atoi(&tmp[j+4]);
			for(j = j+4; j < (int)strlen(tmp); j++) if(tmp[j] == ' ' || tmp[j] == '\t') break;
			tsnp.pos = atoi(&tmp[j+1]);
			tsnp.qual = -1;
			if(STDFORMAT == 1)
			{	for(j = j+1; j < (int)strlen(tmp); j++) if(tmp[j] == ' ' || tmp[j] == '\t') break;
				if(j < (int)strlen(tmp) - 1) tsnp.qual = atof(&tmp[j+1]);
				for(j = j+1; j < (int)strlen(tmp); j++) if(tmp[j] == ' ' || tmp[j] == '\t') break;
			}
			tsnp.alleles.resize(2);
			tsnp.alleles[0] = '0'; tsnp.alleles[1] = '1'; 
			if(STDFORMAT == 1 && j < (int)strlen(tmp) - 3)
			{	tsnp.alleles[0] = tmp[j+1];
				tsnp.alleles[1] = tmp[j+3];
			}
			tsnp.label = 1;
			mydata.snpinfo.push_back(tsnp);
		}
		i++;
	}
	fclose(f);
	mydata.L = (int)mydata.snpinfo.size();

	//read read counts
	double *likP=NULL;
	int c1, c2;
	vector<vector<int> > counts;
	sprintf(fname, "%s.counts", input);
	f = fopen(fname, "r");
	if(f != NULL && !likelihood) 
	{
		if(STDFORMAT)
		{	char tmp[100000], str1[200];
			fgets(tmp, 100000, f);
			int l = (int)strlen(tmp) - 1;
			for(j = 0; j < l; j++) if(tmp[j] == ' ' || tmp[j] == '\t') break;
			for(j = j+1; j < l; j++) if(tmp[j] == ' ' || tmp[j] == '\t') break;
			int N = 0, oj = j + 1;
			do { 	N++;
				for(j = oj; j < l; j++) if(tmp[j] == ' ' || tmp[j] == '\t') break;
				tmp[j] = 0; mydata.indinfo.push_back(&tmp[oj]);
				tmp[j] = ' ';
				oj = j + 1;
			} while(j < l);
			N=N/2;
			counts.resize(N, vector<int>(mydata.L*2, 0));
			for(i = 0; i < mydata.L; i++)
			{	fscanf(f, "%s %s ", &tmp[0], &str1[0]);
				for(j = 0; j < N; j++) fscanf(f, "%d %d ", &counts[j][i*2], &counts[j][i*2+1]);
				fscanf(f, "\n");
			}	
			printf("N=%d ", N);
		}
		else
		{	vector<int> row(mydata.L * 2);
			while(feof(f) == 0)
			{	for(j = 0; j < 5; j++) fscanf(f, "%d ", &k);
				for(j = 0; j < mydata.L; j++) 
				{	fscanf(f, "%d %d ", &c1, &c2);
					row[j * 2] = c1;
					row[j * 2 + 1] = c2;
				}
				counts.push_back(row);
				fscanf(f, "\n");
			}
		}
		fclose(f);
	}
	else if(likelihood)
	{	char tmpname[200];
		sprintf(tmpname, "%s.liks", input);
		f = fopen(tmpname, "r");
		if(f != NULL)
		{
			char tmp[100000], str1[200];
			fgets(tmp, 100000, f);
			int l = (int)strlen(tmp) - 1;
			for(j = 0; j < l; j++) if(tmp[j] == ' ' || tmp[j] == '\t') break;
			for(j = j+1; j < l; j++) if(tmp[j] == ' ' || tmp[j] == '\t') break;
			int N = 0, oj = j + 1;;
			do { 	N++;
				for(j = oj; j < l; j++) if(tmp[j] == ' ' || tmp[j] == '\t') break;
				tmp[j] = 0;
				if(N % 3 != 2) mydata.indinfo.push_back(&tmp[oj]);
				tmp[j] = ' ';
				oj = j + 1;
			} while(j < l);
			N=N/3;
			likP = new double[N * mydata.L * 3];
			char nn1[100], nn2[100], nn3[100];
			for(i = 0; i < mydata.L; i++)
			{	fscanf(f, "%s %s ", &tmp[0], &str1[0]);
				for(j = 0; j < N; j++) 
				{	fscanf(f, "%s %s %s ", &nn1[0], &nn2[0], &nn3[0]);
					*(likP+(j*mydata.L+i)*3) = atof(nn1);
					*(likP+(j*mydata.L+i)*3+1) = atof(nn2);
					*(likP+(j*mydata.L+i)*3+2) = atof(nn3);
				}
				fscanf(f, "\n");
			}	
			printf("N=%d ", N);
			fclose(f);
			counts.resize(N, vector<int>(mydata.L*2, 1));
		}
		else printf("cannot open %s\n", tmpname),exit(0);
	}
	else printf("cannot open %s\n", fname),exit(0);
	mydata.indN = (int)counts.size();//need to multiply ploidity later

	//read fixState
	mydata.fixState.clear();
	if((int)fS.size() > 0) 
	{	if((int)fS.size() != (int)fA.size()) printf("!!!# of haplotype state files does not match with # of haplotype allele files.\n"),exit(0);
		int curK = 0;
		for(int tt = 0; tt < (int)fS.size(); tt++)
		{	vector<vector<int> > tmpfixstate;
			_readPop(fS[tt].c_str(), tmpfixstate);
			int mK = 0;

			for(i = 0; i < (int)tmpfixstate.size(); i++) 
				for(j = 0; j < (int)tmpfixstate[i].size(); j+= 2)
				{	mK = max(mK, tmpfixstate[i][j + 1]);
					tmpfixstate[i][j + 1] += curK;
				}
			mydata.fixState.insert(mydata.fixState.end(), tmpfixstate.begin(), tmpfixstate.end());
			curK += mK + 1;
		}
	}

	//read fixAllele
	mydata.fixAllele.clear();
	vector<vector<char> > falleles;
	int aN = 0;
	if((int)fA.size() > 0)
	{	char tmp[100000];
		FILE *f;
		falleles.resize(mydata.L);
		vector<string> refinfo;	
		for(int tt = 0; tt < (int)fA.size(); tt++)
		{	f = fopen(fA[tt].c_str(), "r");
			if(f == NULL) printf("!!!cannot open file %s\n", fA[tt].c_str()),exit(0);
			int an = 0;
			{ fgets(tmp, 100000, f);
  			  int l = (int)strlen(tmp) - 1;
			  for(j = 0; j < l; j++) if(tmp[j] == ' ' || tmp[j] == '\t') break;
			  for(j = j+1; j < l; j++) if(tmp[j] == ' ' || tmp[j] == '\t') break;
			  for(j = j+1; j < l; j++) if(tmp[j] == ' ' || tmp[j] == '\t') break;
			  int oj = j + 1;
			  do { 	an++;
			 	for(j = oj; j < l; j++) if(tmp[j] == ' ' || tmp[j] == '\t') break;
			 	tmp[j] = 0; refinfo.push_back(&tmp[oj]);
				tmp[j] = ' ';
				oj = j + 1;
			  } while(j < l);
			}
			i = 0;
			int lasti = 0;
			while(fgets(tmp, 100000, f) != NULL)
			{	
				int stl = (int)strlen(tmp);
				for(j = 0; j < stl; j++) if(tmp[j] == ' ' || tmp[j] == '\t') break;
				tmp[j] = 0;
				string rid = tmp;
				tmp[j] = ' ';
		//		for(j = j + 1; j < stl; j++) if(tmp[j] == ' ' || tmp[j] == '\t') break;
				int chr = atoi(&tmp[j + 4]);
				for(j = j + 4; j < stl; j++) if(tmp[j] == ' ' || tmp[j] == '\t') break;
				int pos = atoi(&tmp[j + 1]);
				while(i < mydata.L && (chr > mydata.snpinfo[i].chr || (chr == mydata.snpinfo[i].chr && pos > mydata.snpinfo[i].pos))) i++;
				if(i >= mydata.L) break;
				if(mydata.snpinfo[i].chr != chr || mydata.snpinfo[i].pos != pos) continue;
				for(j = j + 1; j < stl; j++) if(tmp[j] == ' ' || tmp[j] == '\t') break;

				k = 0;
				vector<char> row(an);
				for(j = j + 1; j < stl - 1; j += 2)
				{	if(tmp[j] == mydata.snpinfo[i].alleles[0]) row[k] = 0;
					else if(tmp[j] == mydata.snpinfo[i].alleles[1]) row[k] = 1;
					else row[k] = 2;
					k++;
				}
				if(an > k) row.resize(k);
				an = k;
				for(j = lasti; j < i; j++) falleles[j].resize((int)falleles[j].size() + an, 2);
				falleles[i].insert(falleles[i].end(), row.begin(), row.end());
				lasti = i + 1;
			}
			fclose(f);
			for(j = lasti; j < mydata.L; j++) falleles[j].resize((int)falleles[j].size() + an, 2);
		}
		aN = (int)falleles[0].size() / mydata.ploidity;
		mydata.fixAllele.resize(aN * mydata.ploidity);
		int ll = (int)falleles.size();
		for(j = 0; j < aN * mydata.ploidity; j++)
		{	mydata.fixAllele[j].resize(ll, false);
			for(k = 0; k < ll; k++) mydata.fixAllele[j][k] = (falleles[k][j] == 2);
		}
		if((int)fS.size() > 0 && (int)mydata.fixState.size() != aN * 2) printf("!!!# of individuals with fixed states does not match with # of individuals with fixed alleles.\n"),exit(0);
		mydata.indN += aN;
		mydata.indinfo.insert(mydata.indinfo.begin(), refinfo.begin(), refinfo.end());
	}
	mydata.indIndex = new int[mydata.indN + 1];
	for(i = 0; i < mydata.indN + 1; i++) mydata.indIndex[i] = i;
	mydata.totalN = mydata.indN;

	if((int)mydata.fixState.size() > 0) mydata.fixState.resize(mydata.indN * mydata.ploidity, vector<int>());
	if((int)mydata.fixAllele.size() > 0) mydata.fixAllele.resize(mydata.totalN * mydata.ploidity);

	//calculate Poisson probability
	dSZ = mydata.L;
	dataP = new double[mydata.indN * dSZ * 3];
	nullP = new double[mydata.indN * dSZ * 2];
	errP = new double[dSZ];

	double lambda = 0;
	for(i = 0; i < (int)counts.size(); i++) for(j = 0; j < (int)counts[i].size(); j++) lambda += (double)counts[i][j];
	lambda /= (double)counts.size() * (double)mydata.L;
	double Sden = (double)(mydata.snpinfo[mydata.L - 1].pos - mydata.snpinfo[0].pos) / (double)mydata.L;
	calPoislp(mydata.snpinfo, counts, error, lambda, cut, Sden, dataP, nullP, errP, aN);
	if(likP != NULL)	
	{	double *dpp = dataP + aN * mydata.L * 3, *likpp = likP;
		for(i = aN; i < mydata.indN; i++)
		{	for(j = 0; j < mydata.L * 3; j++, dpp++, likpp++)
			{	*dpp = *likpp;
			//	if(j >= (mydata.L - 1) * 3) printf("%f, ", *likpp);
			}
		}
		for(i = 0; i < mydata.indN * dSZ * 2; i++) nullP[i] = 1e-20;
		for(i = 0; i < dSZ; i++) errP[i] = 1e-20;
		delete likP;
	}
	printf("N = %d(%d), L = %d, lambda=%f\n", mydata.indN, aN, mydata.L, lambda);

	/*
	mydata.indinfo.resize(mydata.totalN);
	for(i = 0; i < mydata.totalN; i++)
	{	char tmp[200];
		if(i < aN) sprintf(tmp, "ref%d", i);
		else sprintf(tmp, "id%d", i - aN);
		mydata.indinfo[i] = tmp;
	}*/

	mydata.data = new float[mydata.totalN * mydata.ploidity * mydata.L];
        for(i = 0; i < aN; i++)
        {       for(j = 0; j < mydata.L; j++)
		{	mydata.data[i * 2 * mydata.L + j] = falleles[j][i * 2];
			mydata.data[(i * 2 + 1) * mydata.L + j] = falleles[j][i * 2 + 1];
		}
        }

	mydata.indN *= mydata.ploidity;//requires preset of ploidity
	return lambda;
}

void hapTensor::calPoislp(vector<SNPINFO> const &snpinfo, vector<vector<int> > const &counts, double error, double lambda, int cut, double Sdensity, double *poisp, double *nullp, double *errp, int aN)
{	double ep = 0;
	int L = (int)counts[0].size() / 2, N = (int)counts.size();
	int i, j, k;
	
	int newcut = 1000000000;
	for(i = 0; i < L; i++)
	{	int n0 = 0, n1 = 0;
		for(j = 0; j < N; j++)
		{	n0 += counts[j][i * 2];
			n1 += counts[j][i * 2 + 1];
		}
		if(newcut > n0) newcut = n0;
		if(newcut > n1) newcut = n1;
	}
	if(newcut > cut) cut = newcut;

	error *= 0.75;
	for(i = 0; i < cut; i++) ep += gsl_ran_poisson_pdf(i, (double)N * lambda * error / 3.);
	ep = (1. - ep) * Sdensity;//1000.; //snp density
	if(error == 0) ep = 0;
       	ep = 1. - pow(1. - ep, 3.);
	ep = ep / (ep + 1.);//proportion of error site
	printf("estimated proportion of FP sites = %f\n", ep);
        
	for(i = 0; i < L; i++)
        {       int n0 = 0, n1 = 0;
                for(j = 0; j < N; j++)
                {       n0 += counts[j][i * 2];
                       	n1 += counts[j][i * 2 + 1];
		}
		double epn = 0;
		for(j = 0; j < cut; j++) epn += gsl_ran_binomial_pdf(j, error / 1., n0 + n1);
		epn = 1. - epn;
                double p0 = gsl_ran_binomial_pdf(n0, ((double)n0 + 0.5) / ((double)n0 + (double)n1 + 1.0), n0 + n1);
                double p1 = gsl_ran_binomial_pdf(min(n0, n1), error / 1., n0 + n1) / (epn + 1e-20);
                double epi = ep * p1 / ((1. - ep) * p0 + ep * p1);
		epi = min(0.0001, epi);
		if(snpinfo[i].qual >= 0) epi = pow(10., -(double)snpinfo[i].qual / 30.);

                double err = ((double)n1 + 0.5) / ((double)(n0 + n1) + 1.0);
		errp[i] = epi;
		double dp0 = gsl_ran_binomial_pdf(n1, error / 1., n0 + n1) / (epn + 1e-20);
	        double dp1 = gsl_ran_binomial_pdf(n0, error / 1., n0 + n1) / (epn + 1e-20);
        	double dp = dp0 / (dp0 + dp1 + 1e-20);
		if(epi > 0.1) printf("%d: epi = %f\terr = %f\tdp = %f | %d %d\n", i, epi, err, dp, n0, n1);
                for(j = 0; j < N; j++)
                {       double sump, tp[3], np[2];
			n0 = counts[j][i * 2];
			n1 = counts[j][i * 2 + 1];
			tp[0] = gsl_ran_binomial_pdf(n0, 1. - max(0.0,error / 1.), n0 + n1);
			tp[1] = gsl_ran_binomial_pdf(n0, 0.5, n0 + n1);
			tp[2] = gsl_ran_binomial_pdf(n0, max(0.0,error / 1.), n0 + n1);
			np[0] = gsl_ran_binomial_pdf(n0, 1. - err, n0 + n1) * dp;
			np[1] = gsl_ran_binomial_pdf(n0, 1. - err, n0 + n1) * (1. - dp);
		
			sump = tp[0] * (1. - epi) + np[0] * epi;
			sump += tp[1] * (1. - epi);
			sump += tp[2] * (1. - epi) + np[1] * epi;	

                        poisp[((j + aN) * L + i) * 3] = tp[0] / sump;
			poisp[((j + aN) * L + i) * 3 + 1] = tp[1] / sump;
			poisp[((j + aN) * L + i) * 3 + 2] = tp[2] / sump;

			nullp[((j + aN) * L + i) * 2] = np[0] / sump;
			nullp[((j + aN) * L + i) * 2 + 1] = np[1] / sump;
		}
		for(j = 0; j < aN; j++)
		{	poisp[(j * L + i) * 3] = poisp[(j * L + i) * 3 + 1] = poisp[(j * L + i) * 3 + 2] = 1./3. * (1. - epi);
			nullp[(j * L + i) * 2] = epi * dp;
			nullp[(j * L + i) * 2 + 1] = epi * (1. - dp);
		}
	}
}

void hapTensor::outputHaplotype(char const *fname, MYDATA const &mydata)
{
	int i, j;
	FILE *f = fopen(fname, "w");
	fprintf(f, "ID CHR POS ");
	for(i = 0; i < (int)mydata.indinfo.size(); i++) fprintf(f, "%s ", mydata.indinfo[i].c_str());
	fprintf(f, "\n");
	for(i = 0; i < mydata.L; i++)
	{	fprintf(f, "%s chr%d %d ", mydata.snpinfo[i].snpid.c_str(), mydata.snpinfo[i].chr, mydata.snpinfo[i].pos);
		for(j = 0; j < mydata.indN; j++) 
			fprintf(f, "%c ", mydata.snpinfo[i].alleles[(int)mydata.data[j * mydata.L + i]]);
		fprintf(f, "\n");
	}
	fclose(f);
}

