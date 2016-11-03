#include "datastructure.h"
#include "genomicTensor.h"

unsigned int rseed = 0;

///////////////////////////////////////////////////
int main(int argc, char *argv[])
{
	int i, j;
	int method = 2, ploidity = 1, minCut = 5, admixN = 30, burnin = 50, mcmc = 50, maxHapK = 0, maxGG = 0, maxPos = 0;
	bool sqc = false, lik = false, samplemaximum = true, outputproportion = false;
	double log2 = -1;
	bool add2 = false, ind = false, indind = false, nb = false;
	double error = 0.01, A = 0, recr = 100., heteroVh = 1., minerr = 0.5, maxerr = 100000000.;
	char const *output = argv[1], *fixPop = NULL, *input = argv[1], *fcov = NULL, *fparam = NULL;
	int startK=0, fixC=0;
	vector<string> fS, fA;
	vector<vector<int> > remove;
	for(i = 2; i < argc; i++)
	{	
		if(strcmp(argv[i], "-rseed") == 0)
		{	rseed = atoi(argv[i + 1]);
			i++;
		}
		else if(strcmp(argv[i], "-o") == 0)
		{	output = argv[i + 1];
			i++;
		}
		else if(strcmp(argv[i], "-P") == 0)
		{	maxPos = atoi(argv[i + 1]);
			i++;
		}
		else if(strcmp(argv[i], "-K") == 0)
		{	maxHapK = atoi(argv[i + 1]);
			i++;
		}
		else if(strcmp(argv[i], "-G") == 0)
		{	maxGG = atoi(argv[i + 1]);
			i++;
		}
		else if(strcmp(argv[i], "-p") == 0)
		{	ploidity = atoi(argv[i + 1]);
			i++;
		}
		else if(strcmp(argv[i], "-c") == 0)
		{	minCut = atoi(argv[i + 1]);
			i++;
		}
		else if(strcmp(argv[i], "-err") == 0)
		{	error = atof(argv[i + 1]);
			i++;
		}
		else if(strcmp(argv[i], "-fixPop") == 0)
		{	fixPop = argv[i + 1];
			i++;
		}
		else if(strcmp(argv[i], "-fixHap") == 0)
		{	if((int)fS.size() != (int)fA.size()) printf("-fixHap option cannot be used together with -fixAllele option.\n"),exit(0);
			fS.push_back(argv[i + 1]);
			fA.push_back(argv[i + 2]);
			i += 2;
		}
		else if(strcmp(argv[i], "-fixAllele") == 0)
		{	if((int)fS.size() > 0) printf("-fixAllele option cannot be used together with -fixHap option.\n"),exit(0);
			fA.push_back(argv[i + 1]);
			i++;
		}
		else if(strcmp(argv[i], "-fixC") == 0)
		{	fixC = atoi(argv[i + 1]);
			i++;
		}
		else if(strcmp(argv[i], "-admixn") == 0)
		{	admixN = atoi(argv[i + 1]);
			i++;
		}
		else if(strcmp(argv[i], "-sqc") == 0)
			sqc = true;
		else if(strcmp(argv[i], "-sample") == 0)
		{	mcmc = atoi(argv[i + 2]);
			burnin = atoi(argv[i + 1]);
			i += 2;
		}
		else if(strcmp(argv[i], "-lik") == 0)
			lik = true;
		else if(strcmp(argv[i], "-A") == 0)
		{	A = atof(argv[i + 1]);
			i ++;
		}
		else if(strcmp(argv[i], "-z") == 0)
		{	fcov = argv[i + 1];
			i++;
		}
		else if(strcmp(argv[i], "-max") == 0)
		{	samplemaximum = true;
		}
		else if(strcmp(argv[i], "-outputprop") == 0)
		{	outputproportion = true;
		}
		else if(strcmp(argv[i], "-rec") == 0)
		{	recr = atof(argv[i + 1]);
			i++;
		}
		else if(strcmp(argv[i], "-log2") == 0) 
		{	log2 = 1.;	
			if(argc > i + 1 && argv[i + 1][0] >= 48 && argv[i + 1][0] < 58) 
			{	log2 = atof(argv[i + 1]);
				i++;
				printf("log2 constant = %f\n", log2);
			}
		}
		else if(strcmp(argv[i], "-minerr") == 0)
		{	minerr = atof(argv[i + 1]);
			i++;
		}
		else if(strcmp(argv[i], "-maxerr") == 0)
		{	maxerr = atof(argv[i + 1]);
			i++;
		}
		else if(strcmp(argv[i], "-add2") == 0) add2 = true;
		else if(strcmp(argv[i], "-param") == 0)
		{	fparam = argv[i + 1];
			i++;
		}
		else if(strcmp(argv[i], "-C") == 0)
		{	startK=atoi(argv[i+1]);
			i++;
		}
		else if(strcmp(argv[i], "-ind") == 0)
		{	ind = true;
		}
		else if(strcmp(argv[i], "-indind") == 0)
		{	indind = true;
		}
		else if(strcmp(argv[i], "-h") == 0)
		{	heteroVh = atof(argv[i + 1]);
			i++;
		}
		else if(strcmp(argv[i], "-remove") == 0)
		{	i++;
			int ind = atoi(argv[i]), ll = (int)strlen(argv[i]) - 1;
			vector<int> row;
			for(j = 1; j < ll; j++) if(argv[i][j] == ':') break;
			do {
				row.push_back(atoi(&argv[i][j + 1]));
				for(j = j + 1; j < ll; j++) if(argv[i][j] == ',') break;
			} while(j < ll);
			if((int)remove.size() <= ind) remove.resize(ind + 1, vector<int>());
			remove[ind] = row;
		}
		else if(strcmp(argv[i], "-nb") == 0)
		{	nb = true;
		}
	}

	time_t tst, ted;
	time(&tst);
	{	genomicTensor gtm(rseed);
		gtm.outputproportion = outputproportion;
		gtm.log2 = log2;
		if(fixC > 0) startK = fixC;
		gtm.run(input, fcov, burnin, mcmc, maxHapK, maxGG, maxPos, A, recr, heteroVh, samplemaximum, output, sqc, fparam, add2, minerr, maxerr, ind, indind, startK, fixC, remove, nb); 
	}
	
	time(&ted);
	printf("total time = %dsec\n", (int)(ted - tst));

	return 0;
}

