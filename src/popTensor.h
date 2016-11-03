#include "datastructure.h"
#include "tensorHMMbase.h"

class popTensor : public tensorHMMbase
{
public:
        popTensor(unsigned int rs=0);
        ~popTensor();

	double *dataP;
	int dSZ, trueploidity;
	vector<vector<int> > trace1, trace2;

	void run(char const *finput, int burnin, int mcmc, int maxHapK, int ploidity, double AA, char const *foutput, char const *fS, int admixN);

protected:
        void _getAsz(MYDATA const &mydata);
        void _updatePara_remove(vector<float> &para, float data, float b, int id, int pos);
        void _updatePara_add(vector<float> &para, float data, float b, int id, int pos);
        void  _sumforwardOne(int id, vector<int> const &breaks, int bst, int bed, TENSORPARA &tpara, int ploidity, float *&dp1, float *&dp2, double const *&rp, int const *&nhp, vector<vector<vector<float> > >::const_iterator &ipm, int const *&aszp, double *&lpp, double *&lmm, double *&lss, int SSZ, int MK, vector<vector<int> > const &ss, vector<vector<int> > const &mlist, double *Vh2);
        int _tracebackOne(int id, vector<int> const &breaks, int bst, int bed, TENSORPARA &tpara, int ploidity, int *&pp, int *&pp2, bool *&rr, bool *&rr2, double const *&rp, int const *&nhp, vector<vector<vector<float> > >::const_iterator &ipm, double *&lpp, double *&lmm, double *&lss, int SSZ, int MK, vector<vector<int> > const &ss, vector<vector<int> > const &mlist, double *Vh2);
	double _getLp(int id, vector<int> const &breaks, int bid, bool flip, int ploidity, TENSORPARA const &tpara, float const *dp, float const *dp2, int const *pp, int const *pp2, vector<vector<vector<float> > >::const_iterator ipm, int const *nhp, int const *aszp);
	double _getLp1(double *dddp, int pop, int state, bool add, bool logscale, vector<vector<vector<float> > >::const_iterator ipm, int const *nhp, int const *aszp, int maxK, double a);
	double _imputeOne(int id, vector<int> const &breaks, int bst, int bed, int flip, int ploidity, TENSORPARA const &tpara, float *dp, float *dp2, int const *pp, int const *pp2, vector<vector<vector<float> > >::const_iterator ipm, int const *nhp, int const *aszp);
	void _updatePrior(vector<bool> const &init, bool oneround = false){};

	void readInput(MYDATA &mydata, char const *input, char const *fS, vector<double> &recomb);
	void _readHaplotype(char const *fname, char *haps, int N, int L);

private:
	int _oneStrandProb(int id, vector<int> const &breaks, int bst, int bed, int strand, TENSORPARA const &tpara, double const *rp, int const *nhp, vector<vector<vector<float> > >::const_iterator ipm, int const *aszp, int MK, double *tP, double *sP, int *path);
};

