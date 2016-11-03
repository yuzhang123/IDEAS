#include <time.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include <math.h>

using namespace std;


int main(int argc, char* argv[])
{

	int i, j, k, l;
	char tmp[1000];
	if(argc > 1)
	{	for(i = 1; i <= 48; i++)
		{	sprintf(tmp, "qdel %d", atoi(argv[1]) + i - 1);
			system(tmp);
		}
	}
	else if(true)//test ideas in different data sets 
	{	int nsz[6] = {3,6,12,24,48,96};
		for(i = 0; i < 6; i++)
		for(j = 1; j < 2; j++)
		for(k = 1; k < 2; k++)
		for(l = 40; l < 43; l+=2) 
		{	FILE *f = fopen("run.script", "w");
			fprintf(f, "#PBS -l nodes=1:ppn=1\n#PBS -q lionxg-yzz2\n#PBS -l pmem=8gb\n#PBS -l walltime=124:00:00\n");
			fprintf(f, "cd $HOME/work/SNPanalysis/source/\n");
			fprintf(f, "Rscript testIDEAS.R %d 10 %d 100 %d %d seqtmp2\n", nsz[i], j, k, l);
			fclose(f);
			sprintf(tmp, "qsub run.script");
			system(tmp);
			time_t s, t = time(&t) + 1;
			do { time(&s); } while(s < t);
		}
	}
	else if(true)//test ideas in different parameters 
	{	int nsz[6] = {2,6,24,96};
		for(i = 0; i < 4; i++)
		for(j = 1; j < 2; j++)
		for(k = 1; k < 2; k++)
		for(l = 0; l < 1; l+=2) 
		{	FILE *f = fopen("run.script", "w");
			fprintf(f, "#PBS -l nodes=1:ppn=1\n#PBS -q lionxg-yzz2\n#PBS -l pmem=8gb\n#PBS -l walltime=124:00:00\n");
			fprintf(f, "cd $HOME/work/SNPanalysis/source/\n");
			fprintf(f, "Rscript sensitivity.R %d 10 %d 100 %d %d seqtmp2\n", nsz[i], j, k, l);
			fclose(f);
			sprintf(tmp, "qsub run.script");
			system(tmp);
			time_t s, t = time(&t) + 1;
			do { time(&s); } while(s < t);
		}
	}

	return 0;
}

