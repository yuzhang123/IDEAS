#include "MixGauss.h"

typedef struct MySortType3 {
        int index;
        double score, weight; 
} MySortType3;

bool operator<(const MySortType3 &a, const MySortType3 &b)
{       return a.score < b.score;
} 

extern unsigned int rseed;

MixGauss::MixGauss()
{
	gsl_rng_env_setup();
        T = gsl_rng_default;
        gammar = gsl_rng_alloc(T);
	if(rseed == 0) gsl_rng_set(gammar, (unsigned int)time(NULL) + 137 * clock());
	else gsl_rng_set(gammar, rseed * 11 + 7);

	lessoneUpdate = true;
	indc = false;
	nb = false;
	gausspara = gaussprior = lambda = NULL;
	minerr = 0.1;
	maxerr = 1000000000.;

	groupcode = coderevmap = NULL;
	xmsz = ymsz = NULL;
	xmap = ymap = NULL;
	mapspace = NULL;
}

MixGauss::~MixGauss(void)
{
	if(gausspara != NULL) delete gausspara;
	if(gaussprior != NULL) delete gaussprior;
	if(lambda != NULL) delete lambda;

	if(groupcode != NULL) delete groupcode;
	if(coderevmap != NULL) delete coderevmap;

	if(xmsz != NULL) delete xmsz;
	if(ymsz != NULL) delete ymsz;

	if(mapspace != NULL) delete mapspace;
	if(xmap != NULL) delete xmap;
	if(ymap != NULL) delete ymap;

	gsl_rng_free(gammar);
}

void MixGauss::addPara(int g, float *yy, float *xx, int id)
{
	int j, k, l;
	double *pp = gausspara + (g * groupn + groupcode[id]) * gausssz;
	pp[0]++;
	int pii = 1, tymsz = ymsz[id], tyvsz = tymsz * (tymsz + 1) / 2, txmsz = xmsz[id], txvsz = txmsz * (txmsz + 1) / 2;
	double *ppp = pp + 1 + tymsz;
	for(j = 0; j < tymsz; j++)
	{	pp[pii + j] += (double)yy[j];
		l = j * (j + 1) / 2;
		for(k = 0; k <= j; k++)
		{	ppp[l + k] += (double)yy[j] * yy[k];
		}
	}
	if(txmsz > 0)
	{	pii = 1 + tymsz + tyvsz;
		ppp = pp + pii + txmsz;
		for(j = 0; j < txmsz; j++)
		{	pp[pii + j] += (double)xx[j];
			l = j * (j + 1) / 2;
			for(k = 0; k <= j; k++)
			{	ppp[l + k] += (double)xx[j] * xx[k];
			}
		}
		ppp += txvsz;
		for(j = 0; j < tymsz; j++)
		{	l = j * txmsz;
			for(k = 0; k < txmsz; k++)
			{	ppp[l + k] += (double)yy[j] * xx[k];
			}
		}
	}
}

void MixGauss::removePara(int g, float *yy, float *xx, int id)
{
	int j, k, l;
	double *pp = gausspara + (g * groupn + groupcode[id]) * gausssz;
	pp[0]--;
	int pii = 1, tymsz = ymsz[id], tyvsz = tymsz * (tymsz + 1) / 2, txmsz = xmsz[id], txvsz = txmsz * (txmsz + 1) / 2;
	double *ppp = pp + 1 + tymsz;
	for(j = 0; j < tymsz; j++)
	{	pp[pii + j] -= (double)yy[j];
		l = j * (j + 1) / 2;
		for(k = 0; k <= j; k++)
		{	ppp[l + k] -= (double)yy[j] * yy[k];
		}
	}
	if(txmsz > 0)
	{	pii = 1 + tymsz + tyvsz;
		ppp = pp + pii + txmsz;
		for(j = 0; j < txmsz; j++)
		{	pp[pii + j] -= (double)xx[j];
			l = j * (j + 1) / 2;
			for(k = 0; k <= j; k++)
			{	ppp[l + k] -= (double)xx[j] * xx[k];
			}
		}
		ppp += txvsz;
		for(j = 0; j < tymsz; j++)
		{	l = j * txmsz;
			for(k = 0; k < txmsz; k++)
			{	ppp[l + k] -= (double)yy[j] * xx[k];
			}
		}
	}
}

void MixGauss::computeLP(float *ydata, float *xdata, int id, double A, double *lp)
{	
	int i, j, k;
	if(lambda == NULL) updateLambda(A);

	double *ll = lambda + lambdasz * groupcode[id], *tlp = lp;
	for(i = 0; i < clustersz; i++, ll += lambdasz * groupn, tlp++)
	{
		(*tlp) = _dmvnorm(ll, ydata, xdata, ymsz[id], xmsz[id]);
	}
}

double MixGauss::_dmvnorm(double *lambdap, float *yp, float *xp, int tymsz, int txmsz)
{
	if(nb)
	{	int i;
		double sum = 0;
		for(i = 0; i < tymsz; i++) sum += yp[i];
		double rt = 0, a = lambdap[tymsz], b = lambdap[tymsz + 1];
		rt += gsl_sf_lngamma(sum + b) + sum * log(a);
		rt -= gsl_sf_lngamma(b) + (sum + b) *  log(1+ a);// + gsl_sf_lngamma(sum + 1);
		for(i = 0; i < tymsz; i++)
		{	rt += yp[i] * log(lambdap[i]+1e-10) - gsl_sf_lngamma(yp[i]+1);
		}
		return(rt);
	}

	int j, k;
	double x[tymsz];
	double *tl = lambdap + 1 + tymsz * (txmsz + 1), m, f = 0, *llp = lambdap + 1, *txp = x;
	for(j = 0; j < tymsz; j++, yp++, txp++)
	{	m = 0;
		(*txp) = (double)(*yp) - (*llp);
		llp++;
		for(k = 0; k < txmsz; k++, llp++) 
		{	(*txp) = (*txp)-(*llp) * (double)xp[k];
		}
		for(k = 0; k <= j; k++, tl++)
		{	m += x[k] * (*tl);
		}
		f += m * m; 
	}
	if(lessoneUpdate) 
	{	f *= (lambdap[lambdasz-1] + (double)tymsz + 3.);
		m = 1. + (*tl) * (*tl);
		tl += txmsz + 1;
		for(j = 0; j < txmsz; j++, tl += txmsz + 1)
		{	double s = *tl;
			for(k = (int)indc * j; k <= j; k++)
				s += (double)xp[k] * tl[k + 1];
			m += s * s;
		}
		f += (double)tymsz * log(m);
	}
	return ((*lambdap)-f)/2.;
}

void MixGauss::updateLambda(double A)
{	int i, j, k;	

	if(lambda != NULL) delete lambda;
	lambda = new double[lambdasz * groupn *  clustersz];
	double *ll = lambda, *rr = gaussprior;
	double *tmpspace = new double[maxymsz * maxymsz * 2 + (maxxmsz+1)*(maxxmsz+1)*2 + (maxxmsz + 1) * maxymsz];
	updateParameterPrior(A);
	for(i = 0;i < clustersz; i++, rr += gausssz)
	{	for(j = 0; j < groupn; j++, ll += lambdasz)
			_getLambdaOne(ll, rr, tmpspace, A, i, coderevmap[j]);
	}
	delete tmpspace;
}

void MixGauss::_getLambdaOne(double *ll, double *rr, double *tmpspace, double priorW, int i, int id)
{	int j, k;	

	double *A = tmpspace; //var y
	double *vL = A + maxymsz * maxymsz; //chol var y
	double *B = vL + maxymsz * maxymsz; //xx
	double *vxL = B + (maxxmsz + 1) * (maxxmsz + 1);
	double *C = vxL + (maxxmsz + 1) * (maxxmsz + 1); //yx
	double *D = C + maxymsz * (maxxmsz + 1); //b'X'X
	double *E = D + maxymsz * (maxxmsz + 1); //b'X'Y+Y'Xb-b'X'Xb
	double *F = E + maxymsz * maxymsz; //
	
	double cc = 1e+10;//1.0/10000.;
	int pii, jj, kk;
	int *tymap = ymap[id], *txmap = xmap[id];
	int tymsz = ymsz[id], txmsz = xmsz[id];
	
	int revxmap[maxxmsz], revymap[maxymsz];
	for(j = 0; j < maxxmsz; j++) revxmap[j] = -1;
	for(j = 0; j < txmsz; j++) revxmap[txmap[j]] = j;
	for(j = 0; j < maxymsz; j++) revymap[j] = -1;
	for(j = 0; j < tymsz; j++) revymap[tymap[j]] = j;


	if(nb)
	{	double sum = 0.1, sum2 = 1;
		for(j = 0; j < tymsz; j++)
		{	ll[j] = rr[1 + tymap[j]] / (rr[1] + 1);
			sum += ll[j];
		}
		int pii = 1 + maxymsz;
		for(j = 0; j < maxyvsz; j++)
		{	int aa = (int)((sqrt(1.+8.*(double)(j+1))-3.)/2.+0.999999);
			int bb = j - aa * (aa + 1) / 2;
			int a = revymap[aa];
			int b = revymap[bb];
			if(a >= 0 && b >= 0)
			{	double f = (*(rr + pii + j)) + (double)(a==b) * priorW;
				sum2 += f;
				if(a!=b) sum2 += f;
			}
		}
		sum2 /= (rr[1] + 1);
if(rr[1] < 1)
{	sum = 0;
	for(j=0;j<tymsz;j++) 
	{	ll[j]=gsl_ran_flat(gammar, 0., 1.); 
		sum += ll[j] / (double)tymsz;
		sum2 += pow(ll[j], 2.) / (double)tymsz;
	}
}
		for(j = 0; j < tymsz; j++) ll[j] /= sum;
		ll[tymsz + 1] = (sum2 - sum * sum - sum) / sum;
		if(ll[tymsz + 1] <= 0) ll[tymsz + 1] = 0.01;
		ll[tymsz] = sum / ll[tymsz + 1];
		if(ll[tymsz] <= 0) ll[tymsz] = 0.01;
printf("%f,%f\n", ll[tymsz], ll[tymsz+1]);fflush(stdout);
			
		return;
	}


	{	//get Y'X
		for(j = 0; j < tymsz; j++)
		{	jj = tymap[j];
			pii = 1 + maxymsz + maxyvsz + maxxmsz + maxxvsz + jj * maxxmsz;
			C[j * (txmsz + 1)] = *(rr + jj + 1);
			for(k = 0; k < txmsz; k++)
			{	kk = txmap[k];
				C[j * (txmsz + 1) + k + 1] = *(rr + pii + kk);
			}
		}

		//get X'X
		pii = 1 + maxymsz + maxyvsz;
		B[0] = (*rr) + priorW + 1e-10;
		for(j = 0; j < txmsz; j++)
		{	jj = txmap[j];
			B[j + 1] = B[(1 + j) * (txmsz + 1)] = *(rr + pii + jj);
		}
		pii += maxxmsz;
		for(j = 0; j < maxxvsz; j++)
		{	int a = revxmap[(int)((sqrt(1.+8.*(double)(j+1))-3.)/2.+0.999999)];
			int b = revxmap[j - a * (a + 1) / 2];
			if(a >= 0 && b >= 0)
			{	B[(a + 1) * (txmsz + 1) + 1 + b] = B[(a + 1) + (txmsz + 1) * (b + 1)] = *(rr + pii + j) + (double)(a==b)/cc;
			}
		}

		vector<double> mrow;
		if((int)modelparameter.size() > i + 1 && (int)modelparameter[i+1].size() > 0) mrow = modelparameter[i + 1];
		else if((int)modelparameter.size() > 0 && (int)modelparameter[0].size() > 0) mrow = modelparameter[0];
		if((int)mrow.size() > 0)
		{	for(j = 0; j < (txmsz + 1); j++)
			{	double td = 0;
				jj = 0; if(j > 0) jj = txmap[j - 1] + 1;
				for(k = 0; k < (txmsz + 1); k++)
				{	kk = 0; if(k > 0) kk = txmap[k - 1] + 1;
					if((int)mrow.size() > kk * 2) td += mrow[kk * 2] * (B[k * (txmsz + 1) + j] - ((double)(kk==0)*1e-10+(double)(kk>0)/cc) * (double)(kk==jj));
				}
				for(k = 0; k < tymsz; k++) D[k * (txmsz + 1) + j] = td;
			}
			for(j = 0; j < tymsz; j++)
			{	double td = 0;
				jj = tymap[j];
				for(k = 0; k < (txmsz + 1); k++)
				{	kk = 0; if(k > 0) kk = txmap[k - 1] + 1;
					if((int)mrow.size() > kk * 2) 
					{	td += mrow[kk * 2] * (C[j * (txmsz + 1) + k] - D[j * (txmsz + 1) + k] / 2.);
					}
				}
				for(k = 0; k < tymsz; k++) F[j * tymsz + k] = td;
			}
			for(j = 0; j < tymsz; j++) 
			{	for(k = 0; k < tymsz; k++) { E[j * tymsz + k] = F[j * tymsz + k] + F[k * tymsz + j]; }
				for(k = 0; k < (txmsz + 1); k++) C[j * (txmsz + 1) + k] -= D[j * (txmsz + 1) + k];
			}

			if((int)mrow.size() > 1) B[0] += mrow[1] - 1e-10;
			for(j = 0; j < maxxvsz; j++)
			{	int a = (int)((sqrt(1.+8.*(double)(j+1))-3.)/2.+0.999999), aa = revxmap[a];
				int b = j - a * (a + 1) / 2, bb = revxmap[b];
				if(aa >= 0 && bb >= 0 && a==b)
				{	double dw = 0;
					if((int)mrow.size() > 3 + a * 2) dw = mrow[3+a*2] - 1./cc;
					B[(aa + 1) * (txmsz + 1) + 1 + bb] += dw;
					B[(aa + 1) + (txmsz + 1) * (bb + 1)] += dw;
				}
			}
		}

		//compute coefficients, step1: ll[1]=(Y'X * chol(X'X^{-1}))
		pii = 1 + tymsz * (txmsz + 1) + (tymsz + 1) * tymsz / 2;
		_cholV(B, txmsz + 1, vxL);
		_invL(vxL, txmsz + 1, ll + pii, false);
		double *lll = ll + 1;
		for(j = 0; j < tymsz; j++)
		{	for(k = 0; k < txmsz + 1; k++, lll++)
			{	*lll = 0.;
				for(int tt = 0; tt <= k; tt++)
				{	(*lll) += C[j * (txmsz + 1) + tt] * (*(ll + pii + k * (txmsz + 1) + tt));	}
			}
		}
			
		//compute coefficients, step2: C=(Y'X * (X'X)^{-1})
		for(j = 0; j < tymsz; j++)
		{	for(k = 0; k < txmsz + 1; k++)
			{	C[j * (txmsz + 1) + k] = 0;
				for(int tt = k; tt < txmsz + 1; tt++)
					C[j * (txmsz + 1) + k] += ll[1 + j * (txmsz + 1) + tt] * ll[pii + tt * (txmsz + 1) + k];
			}
		}
			
		//compute Y variance
		pii = 1 + maxymsz;
		for(j = 0; j < maxyvsz; j++)
		{	int aa = (int)((sqrt(1.+8.*(double)(j+1))-3.)/2.+0.999999);
			int bb = j - aa * (aa + 1) / 2;
			int a = revymap[aa];
			int b = revymap[bb];
			if(a >= 0 && b >= 0)
			{	double f = (*(rr + pii + j)) + (double)(a==b) * priorW;
				if((int)mrow.size() > 0) f -= E[a * tymsz + b];
				for(k = 0; k < txmsz + 1; k++) f -= (*(ll + 1 + a * (txmsz + 1) + k)) * (*(ll + 1 + b * (txmsz + 1) + k));
				A[a * tymsz + b] = A[a + tymsz * b] = f;
				if(indc && a!=b) A[a*tymsz+b]=A[a+tymsz*b]=0;
			}
		}

			
		//store coefficients to the right place
		for(j = 0; j < tymsz * (txmsz + 1); j++)
		{	ll[1 + j] = C[j];
			if((int)mrow.size() > 2 * tymap[(j % (txmsz + 1))]) ll[1 + j] += mrow[tymap[(j % (txmsz + 1))] * 2];
		}
	
		/////////////////
		double s[tymsz];
		for(j = 0; j < tymsz; j++) s[j] = sqrt(A[j*tymsz+j]/(*rr+priorW))/minerr;
		for(j = 0; j < tymsz * tymsz; j++) A[j] /= (min(1.,s[(int)(j/tymsz)]) * min(1.,s[j % tymsz]));
		for(j = 0; j < tymsz; j++) s[j] = sqrt(A[j*tymsz+j]/(*rr+priorW))/maxerr;
		for(j = 0; j < tymsz * tymsz; j++) A[j] /= (max(1.,s[(int)(j/tymsz)]) * max(1.,s[j % tymsz]));

		if(!lessoneUpdate) { for(j=0;j<tymsz*tymsz;j++) A[j]/= ((*rr) + priorW); }

/*		printf("\nn=%d {\n",(int)rr[0]);
		for(j=0;j<tymsz;j++)
		{	for(k=0;k<tymsz;k++) printf("%f,",A[j*tymsz+k]/rr[0]);
			printf("\n");
		}
		printf("}\n");
*/		

		_cholV(A, tymsz, vL);
		_invL(vL, tymsz, A, false);
		double *tt = ll + 1 + tymsz * (txmsz + 1);
		double ldet;

		if(*rr<0) printf("!!!%d:%f ",i,*rr);fflush(stdout);
		if(lessoneUpdate) ldet = -(double)tymsz * log(PI) + (_lgammaM(tymsz, ((*rr) + priorW + (double)tymsz)/2. + 1.5) - _lgammaM(tymsz, ((*rr) + priorW + (double)tymsz)/2. + 1.)) * 2.;
		else ldet = -(double)tymsz * log(2. * PI);
		for(j = 0; j < tymsz; j++)
		{	for(k = 0; k <= j; k++, tt++) (*tt) = A[j * tymsz + k];
			if(A[j * tymsz + j] > 0) 
			{	ldet += log(A[j * tymsz + j]) * 2.;
			}
			else 
			{	ldet += -100000000.;
			}
		}
		ll[0] = ldet;
		ll[lambdasz-1] = *rr+priorW;

		if(ldet < -10000.) 
		{	printf("\n!!!%d:%f, %d\n", i, ldet, (int)gausspara[i*gausssz]);
			for(int xx = 0; xx < lambdasz; xx++) printf("%f, ", ll[xx]);
			printf("\n");
			for(int yy = 0; yy < gausssz; yy++) printf("%f, ", gausspara[i*gausssz+yy]);
			printf("\n");
			for(int zz = 0; zz < gausssz; zz++) printf("%f, ", gaussprior[i*gausssz+zz]);
			printf("\n\n");

		for(j=0;j<tymsz;j++)
		{	for(k=0;k<tymsz;k++) printf("%f,",A[j*tymsz+k]);
			printf("\n");
		}
		for(j=0;j<tymsz;j++) printf("%f,",s[j]);
	
			exit(0);
		}
	}
}

void MixGauss::_cholV(double const *A, int n, double *L)
{	int i, j, k;
	double a, *pl1, *pl2, *pl;
	double const *pa;

	pl = L;
	for(i = 0; i < n * n; i++, pl++) (*pl) = 0;
	for(i = 0; i < n; i++)
	{	pa = A + i * n;
		pl = L + i * n;
		for(j = 0; j <= i; j++, pa++, pl++)
		{	pl1 = L + i * n; pl2 = L + j * n;
			a = 0;
			for(k = 0; k < j; k++) a += (*pl1++) * (*pl2++);
			if(i == j) (*pl) = sqrt(max(0., (*pa) - a));
			else if(L[j * n + j] >= 1e-5) (*pl) = ((*pa) - a) / L[j * n + j];
		}
	}
}

void MixGauss::_invL(double const *L, int n, double *iL, bool transpose)
{	int i, j, k;
	double const *pl;
	double *pil = iL;
	for(i = 0; i < n * n; i++, pil++) (*pil) = 0;
	pl = L; pil = iL;
	for(i = 0; i < n; i++, pl += n + 1, pil += n + 1)
		if((*pl) >= 1e-5) (*pil) = 1. / (*pl);

	for(k = 1; k < n; k++)
		for(i = k; i < n; i++)
		{	double sum = 0;
			pl = L + i * n + i - k;
			pil = iL + (i - k) * n + i - k;
			for(j = i - k; j < i; j++, pl++, pil += n)
				sum += (*pl) * (*pil);
			if((*pl) >= 1e-5) (*pil) = - sum / (*pl);
		}
	if(transpose)
	{	for(i = 0; i < n; i++)
			for(j = i + 1; j < n; j++)
			{	iL[i * n + j] = iL[j * n + i];
				iL[j * n + i] = 0;
			}
	}
}

double MixGauss::_lgammaM(int q, double x)
{
	int i;
	double rt = (double)q * (double)(q - 1) / 4. * log(PI);
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

void MixGauss::rearrangeParameter(int *remap, int newclustersz, float **dataYP, float **dataXP, int totalN, int L, double A)
{
	int i, j, k;
	if(lambda != NULL) { delete lambda; lambda = NULL; }

	double *newgausspara = new double[newclustersz * gausssz * groupn];
	double *newgaussprior = new double[newclustersz * gausssz];
	double msum = 0, tppp = 1. / (1. + A);
	for(i = 0; i < newclustersz; i++)
	{	double *ngp1 = newgausspara + i * gausssz * groupn, *ngp2 = newgaussprior + i * gausssz; 
		double *gp1, *gp2;
		if(i < clustersz)
		{	gp1 = gausspara + remap[i] * gausssz * groupn;
			gp2 = gaussprior + remap[i] * gausssz;

			for(j = 0; j < gausssz * groupn; j++)
			{	*ngp1 = *gp1; ngp1++; gp1++;
			}
			for(j = 0; j < gausssz; j++)
			{	*ngp2 = *gp2; ngp2++; gp2++;
			}
			msum += gp2[0];
		}
		else
		{	for(j = 0; j < gausssz * groupn; j++) ngp1[j] = 0;
			for(j = 0; j < gausssz; j++) ngp2[j] = 0;
			ngp2[0] = max(1e-100, 1. - msum) * (1. - tppp) * pow(tppp, (double)(i-clustersz));

			int kk, ll, pii, dii, ccc=0;
			double maxd = 0;
			do
			{	
				kk = (int)gsl_ran_flat(gammar, 0, (double)totalN);
				ll = (int)gsl_ran_flat(gammar, 0, (double)L);
				dii = ll * ymsz[kk];
				for(j = 0; j < ymsz[kk]; j++) maxd = max(maxd, (double)dataYP[kk][dii + j]);
				ccc++;
			} while(maxd < 2. && i > 0 && ccc < 100000);
			pii = 1;
			for(j = 0; j < ymsz[kk]; j++) ngp2[pii + ymap[kk][j]] = dataYP[kk][dii + j] / 4.;
			pii = 1 + maxymsz + maxyvsz; dii = ll * xmsz[kk];
			for(j = 0; j < xmsz[kk]; j++) ngp2[pii + xmap[kk][j]] = dataXP[kk][dii + j] / 4.;
			
			pii = 1 + maxymsz;
			for(j = 0; j < maxyvsz; j++)
			{	int a = (int)((sqrt(1.+8.*(double)(j+1))-3.)/2.+0.999999);
				int b = j - a * (a + 1) / 2;	
				ngp2[pii + j] = ngp2[1 + a] * ngp2[1 + b] + (double)(a==b)*1./4.;
			}
			pii += maxyvsz + maxxmsz; dii = 1 + maxymsz + maxyvsz;
			for(j = 0; j < maxxvsz; j++)
			{	int a = (int)((sqrt(1.+8.*(double)(j+1))-3.)/2.+0.999999);
				int b = j - a * (a + 1) / 2;	
				ngp2[pii + j] = ngp2[dii + a] * ngp2[dii + b] + (double)(a==b)*1./4.;
			}

			pii += maxxvsz;
			for(j = 0; j < maxymsz * maxxmsz; j++)
				ngp2[pii + j] = ngp2[1 + (int)(j/maxxmsz)] * ngp2[dii + (j % maxxmsz)];
		}
	}
	delete gausspara;
	delete gaussprior;
	gausspara = newgausspara;
	gaussprior = newgaussprior;

	clustersz = newclustersz;
}

void MixGauss::initializePara(int clusterSZ, double *vh, float **dataYP, float **dataXP, int totalN, int L, int maxysz, int *ysz, int **ymp, int maxxsz, int *xsz, int **xmp, double A)
{	int i, j;
	clustersz = clusterSZ;
	maxymsz = maxysz;
	maxyvsz = (maxysz + 1) * maxysz / 2;
	maxxmsz = maxxsz;
	maxxvsz = (maxxsz + 1) * maxxsz / 2;
	maxyxsz = maxysz * maxxsz;
	ymsz = new int[totalN + 1];
	for(i = 0; i < totalN; i++) ymsz[i] = ysz[i];
	xmsz = new int[totalN + 1];
	for(i = 0; i < totalN; i++) xmsz[i] = xsz[i];
	mapspace = new int[(totalN + 1) * (maxymsz + maxxmsz)];
	ymap = new int*[totalN + 1];
	xmap = new int*[totalN + 1];
	for(i = 0; i < totalN + 1; i++) 
	{	ymap[i] = mapspace + i * maxymsz;
		xmap[i] = mapspace + (totalN + 1) * maxymsz + i * maxxmsz;
		if(i < totalN)
		{	for(j = 0; j < ymsz[i]; j++) ymap[i][j] = ymp[i][j];
			for(j = 0; j < xmsz[i]; j++) xmap[i][j] = xmp[i][j];
		}
	}
	_prepareGroupCode(totalN);

	for(i = (int)modelparameter.size() - 1; i > clustersz; i--) if((int)modelparameter[i].size() > 0) { clustersz = i; break; }

	gausssz = 1 + maxymsz + maxyvsz + maxxmsz + maxxvsz + maxyxsz;
	
	lambdasz = 1 + maxymsz * (maxxmsz + 1) + maxyvsz + (maxxmsz + 1) * (maxxmsz + 1) + 1;
	gausspara = new double[clustersz * gausssz * groupn];
	gaussprior = new double[clustersz * gausssz];
	for(i = 0; i < clustersz * gausssz * groupn; i++) gausspara[i] = 0;
	for(i = 0; i < clustersz * gausssz; i++) gaussprior[i] = 0;

	float meany[maxymsz], sdy[maxymsz], meanx[maxxmsz + 1], sdx[maxxmsz + 1], yn[maxymsz], xn[maxxmsz];
	for(i = 0; i < maxymsz; i++) meany[i] = sdy[i] = yn[i] = 0;
	for(i = 0; i < maxxmsz; i++) meanx[i] = sdx[i] = xn[i] = 0;
	for(i = 0; i < totalN; i++)
	{	for(j = 0; j < L * ymsz[i]; j++)
		{	meany[ymap[i][j%ymsz[i]]] += dataYP[i][j];
			sdy[ymap[i][j%ymsz[i]]] += dataYP[i][j] * dataYP[i][j];
		}
		for(j = 0; j < ymsz[i]; j++) yn[ymap[i][j]] += L;
		for(j = 0; j < L * xmsz[i]; j++)
		{	meanx[xmap[i][j%xmsz[i]]] += dataXP[i][j];
			sdx[xmap[i][j%xmsz[i]]] += dataXP[i][j] * dataXP[i][j];
		}
		for(j = 0; j < xmsz[i]; j++) xn[xmap[i][j]] += L;
	}
	for(i = 0; i < maxymsz; i++)
	{	meany[i] /= (yn[i] + 1e-10);
		sdy[i] = sdy[i] / (yn[i] + 1e-10) - meany[i] * meany[i];
	}
	for(i = 0; i < maxxmsz; i++)
	{	meanx[i] /= (xn[i] + 1e-10);
		sdx[i] = sdx[i] / (xn[i] + 1e-10) - meanx[i] * meanx[i];
	}

	for(i = 0; i < clustersz; i++)
	{	gaussprior[i * gausssz] = 1;//vh[i];
		int pii = i * gausssz + 1; 
		for(j = 0; j < maxymsz; j++) gaussprior[pii + j] = meany[j] + gsl_ran_flat(gammar, -sdy[j] * 2., sdy[j] * 2.);
		pii += maxymsz + maxyvsz; 
		for(j = 0; j < maxxmsz; j++) gaussprior[pii + j] = meanx[j] + gsl_ran_flat(gammar, -sdx[j] * 2., sdx[j] * 2.);

		pii -= maxymsz + maxyvsz;
		for(j = 0; j < maxyvsz; j++)
		{	int a = (int)((sqrt(1.+8.*(double)(j+1))-3.)/2.+0.999999);
			int b = j - a * (a + 1) / 2;	
			gaussprior[pii + maxymsz + j] = (gaussprior[pii + a] * gaussprior[pii + b] + (double)(a==b)*1./1. * A);
		}

		pii += maxymsz + maxyvsz;
		for(j = 0; j < maxxvsz; j++) 
		{	int a = (int)((sqrt(1.+8.*(double)(j+1))-3.)/2.+0.999999);
			int b = j - a * (a + 1) / 2;	
			gaussprior[pii + maxxmsz + j] = (gaussprior[pii + a] * gaussprior[pii + b] + (double)(a==b)*1./1. * A);
		}

		pii += maxxmsz + maxxvsz; 
		int opii = i * gausssz + 1;
		for(j = 0; j < maxymsz * maxxmsz; j++)
			gaussprior[pii + j] = gaussprior[opii + (int)(j / maxxmsz)] * gaussprior[opii + maxymsz + maxyvsz + j % maxxmsz];
	}
}

void MixGauss::simData(int n, int id, int g, float *yp)
{	int i, j, k, l;
	double *vl = new double[2 * maxymsz * maxymsz], *ivl = vl + maxymsz * maxymsz;
	double *lambdap = lambda + g * lambdasz * groupn + groupcode[id] * lambdasz, *tl = lambdap + 1 + ymsz[id] * (xmsz[id] + 1);
	for(i = 0; i < ymsz[id]; i++)
	{	for(j = 0; j <= i; j++, tl++)
		{	vl[j * ymsz[id] + i] = 0;
			vl[i * ymsz[id] + j] = (*tl) * (1. + (double)(lessoneUpdate) * (sqrt(lambdap[lambdasz-1]) - 1.));
		}
	}
	_invL(vl, ymsz[id], ivl, false);

	float *typ = yp;
	for(i = 0; i < n; i++)
	{	double z[ymsz[id]], *ttlp = lambdap + 1;;
		for(j = 0; j < ymsz[id]; j++) z[j] = gsl_ran_gaussian(gammar, 1.);
		for(j = 0; j < ymsz[id]; j++, typ++, ttlp += (xmsz[id] + 1))
		{	*typ = (float)*ttlp;
			l = j * ymsz[id];
			for(k = 0; k <= j; k++)
				*typ += (float)(z[k] * ivl[l + k]);
		}
	}
	delete vl;
}

void MixGauss::updateParameterPrior(double A)
{	int i, j;
	double *pp = gausspara, *rr = gaussprior;

	for(i = 0; i < clustersz; i++, pp += gausssz * groupn, rr += gausssz)
	{	_gaussEM(pp, rr, A);
	}
}

void MixGauss::getStateCount(double *cn)
{	int i, j;
	double *pp = gausspara;
	for(i = 0; i < clustersz; i++) 
	{	cn[i] = 0;
		for(j = 0; j < groupn; j++, pp += gausssz)
			cn[i] += *pp;
	}
}

void MixGauss::getStatePrior(double *&ppp, int &step, double A, double B)
{	int i;
	double *rr = gaussprior, tn[clustersz];
	double nn = 0;
	for(i = 0; i < clustersz; i++, rr += gausssz)
	{	tn[i] = /*rr[0];*/max(0., rr[0] - B);
		nn += tn[i];
	}
	_updateVh(tn, clustersz, A, tmpprop);
//	for(i = 0; i < clustersz; i++) tmpprop[i] = tn[i] / nn;
	ppp = &tmpprop[0];
	step = 1;
}

void MixGauss::_updateVh(double *Q, int K, double A, double *V)
{       int i;
	double totalN = 0, *qq = Q, *vv;

	for(i = 0; i < K; i++, qq++) totalN += *qq;

	double vhprod = 1., v;
	qq = Q;
	v = V[0] = (*qq + 1.) / (totalN + 1. + A);
	vv = V + 1;
	for(i = 1; i < K; i++, vv++)
	{       totalN -= *qq;
		vhprod *= 1. - v;
		qq++;
		v = (*qq + 1.) / (totalN + 1. + A);
		*vv = v * vhprod;
	}
	vhprod *= 1. - v;
	V[K] = vhprod;
}

double MixGauss::splitmergeCluster(int type, float **dataYP, float **dataXP, int totalN, int L, float *states, double A, int &mi, int &mj, int &tid)
{	
	mi = mj = tid = -3;
	int i, j, k, l;
	double *gp, cn[clustersz];
	getStateCount(cn);
	int available[clustersz], an = 0;
	for(i = 0; i < clustersz; i++)
	{	if(cn[i] == 0) available[an++] = i;
	}
	if(an == 0) return -10000000000.;
	if(lambda == NULL) updateLambda(A);
	if(type == 0)
	{	double dist[clustersz * (clustersz - 1) / 2], *ll = lambda, *gp1, *ll1;
		gp = gausspara;
		k = 0;
		for(i = 0; i < clustersz - 1; i++, gp += gausssz * groupn, ll += lambdasz * groupn)
		{	gp1 = gp + gausssz * groupn; ll1 = ll + lambdasz * groupn; 
			float my[maxymsz], mx[maxxmsz+1];
			_getMean(gp, &my[0], &mx[0]);

			for(j = i + 1; j < clustersz; j++, gp1 += gausssz * groupn, ll1 += lambdasz * groupn)
			{	double d = 0;
				if(cn[i] <= 0 || cn[j] <= 0) 
				{	d = 100000000.;	}
				else
				{	float my1[maxymsz], mx1[maxxmsz+1], nn = 0;
					_getMean(gp1, &my1[0], &mx1[0]);
					for(int gn = 0; gn < groupn; gn++)
					{	float tmy[maxymsz], tmx[maxxmsz + 1], tmy1[maxymsz], tmx1[maxxmsz + 1];
						int id = coderevmap[gn], kk;
						for(kk = 0; kk < ymsz[id]; kk++)
						{	tmy[kk] = my[ymap[id][kk]];
							tmy1[kk] = my1[ymap[id][kk]];
						}
						for(kk = 0; kk < xmsz[id]; kk++)
						{	tmx[kk] = mx[xmap[id][kk]];
							tmx1[kk] = mx1[xmap[id][kk]];
						}
						d += (min((-_dmvnorm(ll + gn * lambdasz, tmy1, tmx1, ymsz[id], xmsz[id]) + ll[gn * lambdasz] / 2.) / (double)(ymsz[id] + xmsz[id]), (- _dmvnorm(ll1 + gn * lambdasz, tmy, tmx, ymsz[id], xmsz[id]) + ll1[gn * lambdasz] / 2.) / (double)(ymsz[id] + xmsz[id]))) * (gp[gn * gausssz] + gp1[gn * gausssz]); 
						nn += gp[gn * gausssz] + gp1[gn * gausssz];
					}
					d /= nn;
				}	
				dist[k++] = d;// * (min(gp[0], gp1[0]) + 1.);
			}
		}
		k = 0;
		double mind;
		for(i = 0; i < clustersz - 1; i++)
		{	for(j = i + 1; j < clustersz; j++, k++)
			{	if((i == 0 && j == 1) || mind > dist[k])
				{	mind = dist[k]; mi = i; mj = j;	}
//printf("[%d:%d, %d:%d, %f]\n", i, j, (int)gausspara[i*gausssz*groupn], (int)gausspara[j*gausssz*groupn], dist[k]);
			}
		}

		if(mind/* / (double)maxymsz*/ < 1.)
		{	tid = available[0];
			gp = gausspara + mi * gausssz * groupn;
			gp1 = gausspara + mj * gausssz * groupn;
			double *gg = gausspara + tid * gausssz * groupn;
			double *rr = gaussprior + tid * gausssz, *r1 = gaussprior + mi * gausssz, *r2 = gaussprior + mj * gausssz;
			for(i = 0; i < groupn * gausssz; i++, gg++, gp++, gp1++)
				*gg = *gp + *gp1;
		/*	for(i = 0; i < gausssz; i++, rr++, r1++, r2++) 
			{	if(i == 0) *rr = *r1 + *r2;
				else *rr = (*r1 * (double)cn[mi] + *r2 * (double)cn[mj]) / (double)(cn[mi] + cn[mj]);
			}
		*/	
			gg = gausspara + tid * gausssz * groupn;
			gp = gausspara + mi * gausssz * groupn;
			gp1 = gausspara + mj * gausssz * groupn;
			double *lg = lambda + tid * lambdasz * groupn;
			ll = lambda + mi * lambdasz * groupn;
			ll1 = lambda + mj * lambdasz * groupn;
			rr = gaussprior + tid * gausssz;
			double *tmpspace = new double[maxymsz * maxymsz * 4 + (maxxmsz+1)*(maxxmsz+1)*2 + maxymsz * (maxxmsz+1) * 2];
			//updateParameterPrior(A);
			updateLambda(A);


			float *ddd, *xxx;
			float *dp = states;
			double lp1 = 0, lp2 = 0, lp = 0, t1 = 0, t2 = 0;
			for(i = 0; i < totalN; i++)
			{	ddd=dataYP[i];
				xxx=dataXP[i];
				for(j = 0; j < L; j++, dp++, ddd += ymsz[i], xxx += xmsz[i])
				{	if((int)(*dp) == mi)
					{	lp1 += _dmvnorm(ll + groupcode[i] * lambdasz, ddd, xxx, ymsz[i], xmsz[i]);
						double f = _dmvnorm(lg + groupcode[i] * lambdasz, ddd, xxx, ymsz[i], xmsz[i]);
						lp += f;//_dmvnorm(lg + groupcode[i] * lambdasz, ddd, xxx, ymsz[i], xmsz[i]);
						t1 += f;
					}
					else if((int)(*dp) == mj)
					{	lp2 += _dmvnorm(ll1 + groupcode[i] * lambdasz, ddd, xxx, ymsz[i], xmsz[i]);
						double f = _dmvnorm(lg + groupcode[i] * lambdasz, ddd, xxx, ymsz[i], xmsz[i]);
						lp += f;//_dmvnorm(lg + groupcode[i] * lambdasz, ddd, xxx, ymsz[i], xmsz[i]);
						t2 += f;
					}
				}
			}
			delete tmpspace;

//printf("lp=%f, lp1=%f, lp2=%f, t1=%f, t2=%f\n", lp, lp1, lp2, t1, t2);
/*for(i=0;i<gausssz;i++) printf("%f,", rr[i]);printf("\n");
for(i=0;i<gausssz;i++) printf("%f,", r1[i]);printf("\n");
for(i=0;i<lambdasz*groupn;i++) printf("%f,", lg[i]);printf("\n");
for(i=0;i<lambdasz*groupn;i++) printf("%f,", ll[i]);printf("\n");
exit(0);
*/			return lp-lp1-lp2;
		}
		else 
		{	mi = mj = tid = -1;
			return -1000000000.;
		}
	}
	else
	{	
	        if(lambda == NULL) updateLambda(A);
        	gp = gausspara;
 		double sd[clustersz], cumsum[clustersz], md = 0, mdn = 0;
		int sdi[clustersz];
	        for(i = 0; i < clustersz; i++, gp += gausssz * groupn)
	        {      if(cn[i] <= 0) { sd[i] = 0.; sdi[i] = 0; }
                	else
	                {	double tsd[maxymsz], tn[maxymsz];
				for(l = 0; l < maxymsz; l++) tsd[l] = tn[l] = 0;
				for(l = 0; l < groupn; l++)
				{	int tymsz = ymsz[coderevmap[l]], pii = 1 + tymsz;
					k = 0;
					for(j = 0; j < tymsz; j++)
					{	tsd[ymap[coderevmap[l]][j]] += (gp[pii + k] / (gp[0] + 1e-100) - pow(gp[1 + j] / (gp[0] + 1e-100), 2.)) * gp[0];
						tn[ymap[coderevmap[l]][j]] += gp[0];
						k += j + 2;
					}
				}
				for(j = 0; j < maxymsz; j++)
					if(j == 0 || sd[i] < tsd[j] / (tn[j]+1e-100)) { sd[i] = tsd[j] / (tn[j]+1e-100); sdi[i] = j; }
			}	
			md += sd[i] * (double)cn[i];
			mdn += (double)cn[i];
			cumsum[i] = 0;
			if(sd[i]>0) cumsum[i] = exp(sd[i]);
			if(i > 0) cumsum[i] += cumsum[i - 1];
	        }
		md /= mdn;
		double un = gsl_ran_flat(gammar, 0, cumsum[clustersz-1]);
		for(i = 0; i < clustersz; i++) if(un <= cumsum[i]) break;
		mi = i; 
		int mik = sdi[mi];

	        if(sd[mi] > md * 2.)
	        {       tid = available[0];
			double *lg = lambda + mi * lambdasz * groupn, *ll = lambda + tid * lambdasz * groupn;

			float *ddd, *xxx;
        	        float *dp = states;
                	for(i = 0; i < totalN; i++)
	                {       ddd=dataYP[i];
				xxx=dataXP[i];
                        	for(j = 0; j < L; j++, dp++, ddd += ymsz[i], xxx += xmsz[i])
	                        {	if((int)(*dp) == mi)
        	                        {	double ty = ddd[mik] - *(lg + groupcode[i] * lambdasz + 1);
						for(k = 0; k < xmsz[i]; k++) ty -= *(lg + groupcode[i] * lambdasz + 2 + k) * xxx[k];
						if(gsl_ran_flat(gammar, 0., 1.) < exp(ty) / (1. + exp(ty)))//ty > 0) 
						{	*dp = tid;
							removePara(mi, ddd,xxx, i);
							addPara(tid, ddd, xxx, i);
						} 
					}
                	        }
	                }
			
			double oll[lambdasz * groupn];
			for(i = 0; i < lambdasz * groupn; i++) oll[i] = lg[i];
                	double *tmpspace = new double[maxymsz * maxymsz * 4 + (maxxmsz+1)*(maxxmsz+1)*2 + maxymsz * (maxxmsz+1) * 2];
			double *gp = gausspara + mi * gausssz * groupn, *gg = gausspara + tid * gausssz * groupn;
                	double *r1 = gaussprior + mi * gausssz, *rr = gaussprior + tid * gausssz;
			//updateParameterPrior(A);
			updateLambda(A);

        	        dp = states;
	                double lp1 = 0, lp2 = 0, lp = 0;
        	        for(i = 0; i < totalN; i++)
                	{       ddd=dataYP[i];
				xxx=dataXP[i];
	                        for(j = 0; j < L; j++, dp++, ddd += ymsz[i], xxx += xmsz[i])
        	                {       if((int)(*dp) == mi)
                	                {       lp1 += _dmvnorm(lg + groupcode[i] * lambdasz, ddd, xxx, ymsz[i], xmsz[i]);
                        	                lp += _dmvnorm(&oll[0] + groupcode[i] * lambdasz, ddd, xxx, ymsz[i], xmsz[i]);
                                	}
	                                else if((int)(*dp) == tid)
        	                        {       lp2 += _dmvnorm(ll + groupcode[i] * lambdasz, ddd, xxx, ymsz[i], xmsz[i]);
                	                        lp += _dmvnorm(&oll[0] + groupcode[i] * lambdasz, ddd, xxx, ymsz[i], xmsz[i]);
                        	        }
                        	}
			}
			delete tmpspace;
			return lp-lp1-lp2;
		}
		else 
		{	mi = mj = tid = -1;
			return -1000000000.;
		}
	}
}

void MixGauss::clearParameter()
{	int i;
	if(gausspara != NULL) delete gausspara;
	if(gaussprior != NULL) delete gaussprior;
	gausspara = new double[clustersz * gausssz * groupn];
	gaussprior = new double[clustersz * gausssz];
	for(i = 0; i < clustersz * gausssz * groupn; i++) gausspara[i] = 0;
	for(i = 0; i < clustersz * gausssz; i++) gaussprior[i] = 0;
	if(lambda != NULL)
	{	delete lambda;
		lambda = NULL;
	}
}

void MixGauss::outputParameter(char *fout, double A)
{	int i, j;
	FILE *f = fopen(fout, "w");
	for(i = 0; i < clustersz; i++)
	{
		fprintf(f, "%f\t", gaussprior[i * gausssz]);
		for(j = 1; j < gausssz; j++)
			fprintf(f, "%f\t", gaussprior[i * gausssz + j]);
		fprintf(f, "\n");
	}
	fclose(f);
}

void MixGauss::_prepareGroupCode(int totalN)
{	if(groupcode != NULL) delete groupcode;
	if(coderevmap != NULL) delete coderevmap;

	groupcode = new int[totalN + 1];
	coderevmap = new int[totalN + 1];
	int i, j, k, l;
	groupn = 0;
	for(i = 0; i < totalN; i++)
	{	for(j = 0; j < groupn; j++)
		{	if(ymsz[i] == ymsz[coderevmap[j]] && xmsz[i] == xmsz[coderevmap[j]])
			{	for(k = 0; k < ymsz[i]; k++)	
					if(ymap[i][k] != ymap[coderevmap[j]][k]) break;
				if(k >= ymsz[i]) 
				{	for(l = 0; l < xmsz[i]; l++)
						if(xmap[i][l] != xmap[coderevmap[j]][l]) break;
					if(l >= xmsz[i]) break;
				}
			}
		}
		if(j >= groupn)
		{	coderevmap[groupn] = i;
			groupn++;
		}
		groupcode[i] = j;
	}
/*
printf("groupn=%d ", groupn);
for(i=0;i<totalN;i++) printf("%d,",groupcode[i]);
for(i=0;i<groupn;i++) printf("%d,",coderevmap[i]);
printf("\n");
*/
	lambdasz = 1 + maxymsz * (maxxmsz + 1) + maxyvsz + (maxxmsz + 1) * (maxxmsz + 1);
}

void MixGauss::_getMean(double const *pp, float *my, float *mx)
{	int j, l;
	float dy[maxymsz], ny[maxymsz], dx[maxxmsz + 1], nx[maxxmsz + 1];
	for(j = 0; j < maxymsz; j++) dy[j] = ny[j] = 0;
	for(j = 0; j < maxxmsz; j++) dx[j] = nx[j] = 0;
	for(j = 0; j < groupn; j++) 
	{	int ii = coderevmap[j];
		for(l = 0; l < ymsz[ii]; l++)
		{	dy[ymap[ii][l]] += pp[1 + l];
			ny[ymap[ii][l]] += pp[0];
		}
	}
	for(j = 0; j < maxymsz; j++) my[j] = dy[j] / (ny[j] + 1e-100);

	for(j = 0; j < groupn; j++) 
	{	int ii = coderevmap[j], pii = 1 + ymsz[ii] + ymsz[ii] * (ymsz[ii] + 1) / 2;
		for(l = 0; l < xmsz[ii]; l++)
		{	dx[xmap[ii][l]] += pp[pii + l];
			nx[xmap[ii][l]] += pp[0];
		}
	}
	for(j = 0; j < maxxmsz; j++) mx[j] = dx[j] / (nx[j] + 1e-100);
}

void MixGauss::_gaussEM(double *pp, double *rr, double AA)
{	int i, j, k, len = maxymsz + maxxmsz + 1;

	double *opp = pp;
	double *tmpspace = new double[len * len * 7];
	double *phi = tmpspace, *newphi = phi + len * len, *T = newphi + len * len, *tmpuse = T + len * len, datasz; 
	for(i = 0; i < len * len; i++) phi[i] = (double)((int)(i / len) == (i % len)) * AA;

	bool a_inner;
	int round = 0;
	double d, tr;
	do {
		a_inner = true;
		for(i = 0; i < len * len; i++) newphi[i] = (double)((int)(i/len)==(i%len)) * AA;
		pp = opp;

		//compute 2nd moments by EM
		for(i = 0; i < groupn; i++, pp += gausssz)
		{
			if(pp[0] <= 0.5) continue;
			int id = coderevmap[i];
			int tymsz = ymsz[id], tvmsz = tymsz * (tymsz + 1) / 2, txmsz = xmsz[id], txvsz = txmsz * (txmsz + 1) / 2;;

			int misn, obsn, misyn, misxn, misy[maxymsz], misx[maxxmsz + 1], map[len];
			_collapsedPrediction(pp, phi, id, T, tmpuse, misyn, misy, misxn, misx, map);
			misn = misyn + misxn; obsn = tymsz + txmsz + 1;
			double *A = T, *B = A + misn * misn, *C = B + obsn * misn;

			//get var	
			for(j = 0; j < len * (len + 1) / 2; j++)
			{	double v = 0;
				int a = (int)((sqrt(1.+8.*(double)(j+1))-3.)/2.+0.999999);
				int b = j - a * (a + 1) / 2;
				int aa, bb;
				if(a < misn)//mis/mis
					v = A[a * misn + b];
				else if(b < misn)//mis/obs
					v = B[(a - misn) * misn + b];
				else //obs/obs
					v = C[(a - misn) * obsn + (b - misn)];
				aa = map[a]; bb = map[b];
				newphi[aa * len + bb] += v;
				if(aa != bb) newphi[bb * len + aa] += v;
/*if(aa==len-1 && bb==4) 
{	int ll;
	for(ll = 0; ll < misn; ll++) if(map[ll] == 4) break;
	double o=0;
	if(ll>=misn) { for(ll = 0; ll < ymsz[id]; ll++) if(ymap[id][ll]==4) break; if(ll>=ymsz[id]) printf("!!!"),exit(0); o=pp[1+ll]; }
	printf("[%d:%f:%f/%f|", i, o, v, C[obsn * obsn - 1]);
}*/
			}
		}
//printf("\n");

		datasz = newphi[len * len - 1];
		for(i = 0; i < len * len; i++) newphi[i] /= datasz;	
		//check convergence
		d = tr = 0;
		for(i = 0; i < len - 1; i++) 
		{	d += fabs(newphi[i * len + i] - phi[i * len + i]);
			tr += phi[i * len + i];
		}
		//d = sqrt(d / (double)(len - 1));
		//tr /= (double)(len - 1);
		if(d / tr < 0.01) a_inner = false;
/*
printf("<%d: %f/%f=%f>\n", round, d, tr, d/tr);

		for(i = 0; i < len; i++)
		{	for(j = 0; j < len; j++)
				printf("%f ", newphi[i * len + j]);
			printf("\n");
		}; printf("\n");
*/
		for(i = 0; i < len * len; i++) phi[i] = newphi[i];
		round++;
	} while(a_inner && round < 100);
	if(round >= 100) printf("\n!!!warning: EM not converged, d/tr=%f/%f=%f.\n",d,tr,d/tr),fflush(stdout);
/*
pp=opp;
if(phi[4*len+len-1]>2.)
{	printf("\n<%f: ", phi[4*len+len-1]);
	for(i=0;i<groupn;i++)
	{	for(j=0;j<ymsz[coderevmap[i]];j++)
			if(ymap[coderevmap[i]][j] == 4) break;
		if(j < ymsz[coderevmap[i]])
		{	printf("%d:%f,", i,pp[i*gausssz+1+j]);//pp[i*gausssz]);
		}
	}
	printf(">\n");
	k=(int)(opp-gausspara)/gausssz/groupn;
	printf("k=%d\n", k);
	for(i=0;i<43;i++) if(groupcode[i]==3) printf("%d,",i);
exit(0);
}*/
	//fill in rr
	rr[0] = datasz;
//	if(datasz > AA)
	{	int pii = 1;
		for(i = 0; i < maxymsz; i++) rr[pii + i] = phi[i * len + len - 1] * datasz;
		pii = 1 + maxymsz + maxyvsz;
		for(i = 0; i < maxxmsz; i++) rr[pii + i] = phi[(i + maxymsz) * len + len - 1] * datasz;
		pii = 1 + maxymsz;
		for(i = 0; i < maxyvsz; i++)
		{	int a = (int)((sqrt(1.+8.*(double)(i+1))-3.)/2.+0.999999);
			int b = i - a * (a + 1) / 2;
			rr[pii + i] = phi[a * len + b] * datasz;
		}
		pii = 1 + maxymsz + maxyvsz + maxxmsz;
		for(i = 0; i < maxxmsz; i++)
		{	int a = (int)((sqrt(1.+8.*(double)(i+1))-3.)/2.+0.999999);
			int b = i - a * (a + 1) / 2;
			rr[pii + i] = phi[(a + maxymsz) * len + b + maxymsz] * datasz;
		}
		pii = 1 + maxymsz + maxyvsz + maxxmsz + maxxvsz;
		for(i = 0; i < maxyxsz; i++)
		{	int a = (int)(i / maxxmsz);
			int b = i % maxxmsz;
			rr[pii + i] = phi[(b + maxymsz) * len + a] * datasz;
		}
	}
//for(i=0;i<gausssz;i++) printf("%f,",rr[i]);printf("\n");

	delete tmpspace;
}

void MixGauss::_collapsedPrediction(double const *pp, double const *phi, int id, double *rt, double *tmpspace, int &misyn, int *misy, int &misxn, int *misx, int *map)
{
	int i, j, k;
	int tymsz = ymsz[id], tyvsz = tymsz * (tymsz + 1) / 2, txmsz = xmsz[id], txvsz = txmsz * (txmsz + 1) / 2;;
	misyn = maxymsz - tymsz; misxn = maxxmsz - txmsz; 
	int misn = misyn + misxn, obsn = tymsz + txmsz + 1;
	for(j = 0; j < maxymsz; j++) misy[j] = 0;
	for(j = 0; j < tymsz; j++) misy[ymap[id][j]] = 1;
	k = 0; for(j = 0; j < maxymsz; j++) if(misy[j] == 0) misy[k++] = j;
	for(j = 0; j < maxxmsz; j++) misx[j] = 0;
	for(j = 0; j < txmsz; j++) misx[xmap[id][j]] = 1;
	k = 0; for(j = 0; j < maxxmsz; j++) if(misx[j] == 0) misx[k++] = j;

	int len = maxymsz + maxxmsz + 1;
	for(i = 0; i < misyn; i++) map[i] = misy[i];
	for(i = 0; i < misxn; i++) map[misyn + i] = misx[i] + maxymsz;
	for(i = 0; i < tymsz; i++) map[misn + i] = ymap[id][i];
	for(i = 0; i < txmsz; i++) map[misn + tymsz + i] = xmap[id][i] + maxymsz;
	map[len - 1] = len - 1;

	double *A = rt, *B = A + misn * misn, *C = B + obsn * misn;
	double *OO = tmpspace, *MO = OO + obsn * obsn, *MM = MO + obsn * misn, *Beta = MM + misn * misn;
	double *vl = Beta + obsn * misn, *ivl = vl + obsn * obsn;
	for(i = misn; i < len; i++)
		for(j = 0; j < len; j++)
		{	if(i >= misn && j >= misn) OO[(i - misn) * obsn + (j - misn)] = phi[map[i] * len + map[j]];
			else if(j < misn) MO[(i - misn) * misn + j] = phi[map[i] * len + map[j]];
		}
	
	//compute beta=Y'X(X'X)^-1
	_cholV(OO, obsn, vl);
	_invL(vl, obsn, ivl, false);
	_XtX(ivl, obsn, obsn, OO); 
	_AtB(MO, OO, misn, obsn, obsn, Beta);//Beta: misn x obsn matrix
	_AB(Beta, MO, misn, misn, obsn, MM);
/*
if(misn > 0)
{
printf("\nbeta(%d,%d):\n",(int)((pp-gausspara)/gausssz/groupn),(int)(((pp-gausspara)/gausssz)%groupn));
for(i=0;i<misn;i++)
{	for(j=0;j<obsn;j++) printf("%f,",Beta[i*obsn+j]);
	printf("\n");
};printf("\n");
}*/

	//put observations into C
	for(i = 0; i < obsn; i++)
		for(j = 0; j <= i; j++)
		{	double v = 0;
			if(i < tymsz) v = pp[1 + tymsz + i * (i + 1) / 2 + j];
			else if(j >= tymsz)//xx
			{	int ii = i - tymsz, jj = j - tymsz;
				if(i - tymsz < txmsz) v = pp[1 + tymsz + tyvsz + txmsz + ii * (ii + 1) / 2 + jj];
				else if(j - tymsz < txmsz) v = pp[1 + tymsz + tyvsz + jj];
				else v = pp[0];
			}
			else//yx
			{	int ii = i - tymsz;
				if(ii < txmsz) v = pp[1 + tymsz + tyvsz + txmsz + txvsz + j * txmsz + ii];
				else v = pp[1 + j];
			}
			C[i * obsn + j] = C[j * obsn + i] = v;
		}
	_ABt(C, Beta, obsn, misn, obsn, B);
	_AB(Beta, B, misn, misn, obsn, A);
if(A[0]<-0.1)
{	printf("\nmisn=%d, obsn=%d \n", misn, obsn);
	for(i=0;i<gausssz;i++) printf("%f,",pp[i]);printf("\n");
	printf("C:\n");
	for(i=0;i<obsn;i++)
	{	for(j=0;j<obsn;j++) printf("%f,",C[i*obsn+j]);
		printf("\n");
	}
	printf("beta:\n");
	for(i=0;i<misn;i++)
	{	for(j=0;j<obsn;j++) printf("%f,",Beta[i*obsn+j]);
		printf("\n");
	}
	printf("B:\n");
	for(i=0;i<obsn;i++)
	{	for(j=0;j<misn;j++) printf("%f,",B[i*misn+j]);
		printf("\n");
	}
	printf("A:\n");
	for(i=0;i<misn;i++)
	{	for(j=0;j<misn;j++) printf("%f,",A[i*misn+j]);
		printf("\n");
	}
	
	double n = 0;
	printf("gausspara:\n");
	for(i=0;i<clustersz;i++)
	{	printf("i=%d:\n", i);
		for(j = 0; j < groupn; j++)
		{	for(int k = 0; k < gausssz; k++) printf("%f,", gausspara[i * gausssz * groupn + j * gausssz + k]); printf("\n");
			n += gausspara[(i * groupn + j) * gausssz];
		}
	}
	printf("n=%f\n", n);
	exit(0);
}
	for(i = 0; i < misn; i++)
		for(j = 0; j < misn; j++)
			A[i * misn + j] += (phi[map[i] * len + map[j]] - MM[i * misn + j]) * pp[0];
}

void MixGauss::_XtX(double const *X, int rn, int cn, double *X2)
{       int i, j, k;

        double const *pi, *pj;
        for(i = 0; i < cn; i++)
        {       for(j = 0; j <= i; j++)
                {       pi = X + i;
                        pj = X + j;
                        double s = 0;
                        for(k = 0; k < rn; k++, pi += cn, pj += cn)
                                s += (*pi) * (*pj);
                        X2[i * cn + j] = X2[j * cn + i] = s;
                }
        }
}

void MixGauss::_AtB(double const *A, double const *B, int ca, int cb, int n, double *C)
{	int i, j, k;
	double const *pa, *pb;
	double *pc = C;
	for(i = 0; i < ca; i++)
		for(j = 0; j < cb; j++, pc++)
		{	pa = A + i;
			pb = B + j;
			(*pc) = 0;
			for(k = 0; k < n; k++, pa+=ca, pb+=cb)
				(*pc) += (*pa) * (*pb);
		}
}

void MixGauss::_ABt(double const *A, double const *B, int ra, int rb, int n, double *C)
{	int i, j, k;
	double const *pa, *pb;
	double *pc = C;
	for(i = 0; i < ra; i++)
		for(j = 0; j < rb; j++, pc++)
		{	pa = A + i * n;
			pb = B + j * n;
			(*pc) = 0;
			for(k = 0; k < n; k++, pa++, pb++)
				(*pc) += (*pa) * (*pb);
		}
}

void MixGauss::_AB(double const *A, double const *B, int ra, int cb, int n, double *C)
{	int i, j, k;
	double const *pa, *pb;
	double *pc = C;
	for(i = 0; i < ra; i++)
		for(j = 0; j < cb; j++, pc++)
		{	pa = A + i * n;
			pb = B + j;
			(*pc) = 0;
			for(k = 0; k < n; k++, pa++, pb += cb)
				(*pc) += (*pa) * (*pb);
		}
}

void MixGauss::imputeMissingTrack(int id, int L, float *yp, float *xp, float *yimp, float *ximp, float *states)//double *stateprob)
{
	int i, j, k;
	int tymsz = ymsz[id], tyvsz = tymsz * (tymsz + 1) / 2, txmsz = xmsz[id], txvsz = txmsz * (txmsz + 1) / 2;;
	int misyn = maxymsz - tymsz, misxn = maxxmsz - txmsz, misy[maxymsz], misx[maxxmsz + 1]; 
	int misn = misyn + misxn, obsn = tymsz + txmsz + 1, len = maxymsz + maxxmsz + 1, mapy[maxymsz], mapx[maxxmsz + 1];
	if(misn == 0) return;
	for(j = 0; j < maxymsz; j++) misy[j] = 0;
	for(j = 0; j < tymsz; j++) misy[ymap[id][j]] = 1;
	k = 0; for(j = 0; j < maxymsz; j++) if(misy[j] == 0) misy[k++] = j;
	for(j = 0; j < maxxmsz; j++) misx[j] = 0;
	for(j = 0; j < txmsz; j++) misx[xmap[id][j]] = 1;
	k = 0; for(j = 0; j < maxxmsz; j++) if(misx[j] == 0) misx[k++] = j;
	for(i = 0; i < misyn; i++) mapy[misy[i]] = i;
	for(i = 0; i < misxn; i++) mapx[misx[i]] = misyn + i;
	for(i = 0; i < tymsz; i++) mapy[ymap[id][i]] = misn + i;
	for(i = 0; i < txmsz; i++) mapx[xmap[id][i]] = misn + tymsz + i;

	double *tmpspace = new double[obsn * obsn * 3 + obsn * misn * (clustersz + 1)];
	double *OO = tmpspace, *MO = OO + obsn * obsn, *Beta = MO + obsn * misn;
	double *vl = Beta + obsn * misn * clustersz, *ivl = vl + obsn * obsn;
	double *rr = gaussprior;
	for(k = 0; k < clustersz; k++, rr += gausssz)
	{	int a, b, pii;
		pii = 1;
		for(i = 0; i < maxymsz; i++)
		{	a = mapy[i];
			if(a < misn) MO[(obsn-1)*misn + a] = rr[pii + i];
			else OO[(obsn-1)*obsn+a-misn] = OO[(a-misn)*obsn+obsn-1] = rr[pii + i];
		}
		pii = 1 + maxymsz + maxyvsz;
		for(i = 0; i < maxxmsz; i++)
		{	a = mapx[i];
			if(a < misn) MO[(obsn-1)*misn + a] = rr[pii + i];
			else OO[(obsn-1)*obsn+a-misn] = OO[(a-misn)*obsn+obsn-1] = rr[pii + i];
		}
		pii = 1 + maxymsz;
		for(i = 0; i < maxyvsz; i++)
		{	a = (int)((sqrt(1.+8.*(double)(i+1))-3.)/2.+0.999999);
			b = i - a * (a + 1) / 2;
			a = mapy[a]; b = mapy[b];
			if(a < misn && b < misn);
			else if(a >= misn && b >= misn) OO[(a-misn)*obsn+b-misn] = OO[(b-misn)*obsn + a-misn] = rr[pii + i];
			else if(a >= misn && b < misn) MO[(a-misn)*misn + b] = rr[pii + i];
			else MO[(b-misn)*misn+a] = rr[pii + i];
		}
		pii = 1 + maxymsz + maxyvsz + maxxmsz;
		for(i = 0; i < maxxvsz; i++)
		{	a = (int)((sqrt(1.+8.*(double)(i+1))-3.)/2.+0.999999);
			b = i - a * (a + 1) / 2;
			a = mapx[a]; b = mapx[b];
			if(a < misn && b < misn);
			else if(a >= misn && b >= misn) OO[(a-misn)*obsn+b-misn] = OO[(b-misn)*obsn + a-misn] = rr[pii + i];
			else if(a >= misn && b < misn) MO[(a-misn)*misn + b] = rr[pii + i];
			else MO[(b-misn)*misn + a] = rr[pii + i];
		}
		pii = 1 + maxymsz + maxyvsz + maxxmsz + maxxvsz;
		for(i = 0; i < maxymsz; i++)
			for(j = 0; j < maxxmsz; j++)
			{	a = mapy[i]; b = mapx[j];
				if(a < misn && b < misn);
				else if(a >= misn && b >= misn) OO[(a-misn)*obsn+b-misn] = OO[(b-misn)*obsn + a-misn] = rr[pii + i];
				else if(a < misn) MO[(b-misn)*misn+a] = rr[pii + i];
				else MO[(a-misn)*misn+b] = rr[pii + i];
			}
		OO[(obsn-1)*obsn+obsn-1] = rr[0];
	
		//compute beta=Y'X(X'X)^-1
		_cholV(OO, obsn, vl);
		_invL(vl, obsn, ivl, false);
		_XtX(ivl, obsn, obsn, OO); 
		_AtB(MO, OO, misn, obsn, obsn, Beta + misn * obsn * k);//Beta: misn x obsn matrix

if(id==5)
{
		printf("\n%d:\n",k);
/*		printf("MO:\n");
		for(i=0;i<obsn;i++)
		{	for(j=0;j<misn;j++) printf("%f ", MO[i*misn+j]);
			printf("\n");
		}
		printf("OO:\n");
		for(i=0;i<obsn;i++)
		{	for(j=0;j<obsn;j++) printf("%f ", OO[i*obsn+j]);
			printf("\n");
		}
*/		printf("Beta:\n");
		for(i=0;i<misn;i++)
		{	for(j=0;j<obsn;j++) printf("%f ", Beta[misn*obsn*k+i*obsn+j]);
			printf("\n");
		}
		printf("rr: ");
		for(i=0;i<gausssz;i++) printf("%f,",rr[i]);printf("\n");
		for(i=0;i<maxymsz;i++) printf("%d,",mapy[i]);printf("\n");
}
	}

	double pred[obsn], resp[misn], wresp[misn], wt;
	pred[obsn - 1] = 1.;
	for(i = 0; i < L; i++, yp += tymsz, xp += txmsz, yimp += misyn, ximp += misxn, states++)//stateprob += clustersz)
	{	for(j = 0; j < tymsz; j++) pred[j] = yp[j];
		for(j = 0; j < txmsz; j++) pred[tymsz + j] = xp[j];
		
		_AB(Beta + misn * obsn * (int)(*states), &pred[0], misn, 1, obsn, &resp[0]);
		for(j = 0; j < misyn; j++) yimp[j] = (float)resp[j];//,printf("%f:%f ",*states,resp[j]);
		for(j = 0; j < misxn; j++) ximp[j] = (float)resp[misyn + j];
	/*	
		wt = 1e-100;
		for(j = 0; j < misn; j++) wresp[j] = 0;
		for(j = 0; j < clustersz; j++)
		{	_AB(Beta + misn * obsn * j, &pred[0], misn, 1, obsn, &resp[0]);
			for(k = 0; k < misn; k++) wresp[j] += resp[j] * stateprob[j];
			wt += stateprob[j];
		}
		for(j = 0; j < misyn; j++) yimp[j] = (float)wresp[j] / wt;//,printf("%f:%f ",*states,resp[j]);
		for(j = 0; j < misxn; j++) ximp[j] = (float)wresp[misyn + j] / wt;
	*/	
//		for(j=0;j<obsn;j++) printf("%f,", pred[j]); printf("\t");
//		for(j=0;j<misn;j++) printf("%f,", resp[j]); printf("\t%d\n",(int)(*states));
	}

	delete tmpspace;
}
