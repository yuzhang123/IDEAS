#include "datastructure.h"
#include "tensorHMMbase.h"

class hapTensor : public tensorHMMbase
{
public:
        hapTensor(unsigned int rs=0);
        ~hapTensor();

	int trueploidity, dSZ;
	double *dataP, *nullP, *errP;

	void run(char const *finput, bool likelihood, int burnin, int mcmc, int maxHapK, int ploidity, double err, double AA, bool samplemax, int cut, char const *foutput, vector<string> const &fS, vector<string> const &fA, bool sqc);

protected:
        void _getAsz(MYDATA const &mydata);
        void _updatePara_remove(vector<float> &para, float data, float b, int id, int pos);
        void _updatePara_add(vector<float> &para, float data, float b, int id, int pos);
        void  _sumforwardOne(int id, vector<int> const &breaks, int bst, int bed, TENSORPARA &tpara, int ploidity, float *&dp1, float *&dp2, double const *&rp, int const *&nhp, vector<vector<vector<float> > >::const_iterator &ipm, int const *&aszp, double *&lpp, double *&lmm, double *&lss, int SSZ, int MK, vector<vector<int> > const &ss, vector<vector<int> > const &mlist, double *Vh2);
        int _tracebackOne(int id, vector<int> const &breaks, int bst, int bed, TENSORPARA &tpara, int ploidity, int *&pp, int *&pp2, bool *&rr, bool *&rr2, double const *&rp, int const *&nhp, vector<vector<vector<float> > >::const_iterator &ipm, double *&lpp, double *&lmm, double *&lss, int SSZ, int MK, vector<vector<int> > const &ss, vector<vector<int> > const &mlist, double *Vh2);
	double _getLp(int id, vector<int> const &breaks, int bid, bool flip, int ploidity, TENSORPARA const &tpara, float const *dp, float const *dp2, int const *pp, int const *pp2, vector<vector<vector<float> > >::const_iterator ipm, int const *nhp, int const *aszp);
	double _imputeOne(int id, vector<int> const &breaks, int bst, int bed, int flip, int ploidity, TENSORPARA const &tpara, float *dp, float *dp2, int const *pp, int const *pp2, vector<vector<vector<float> > >::const_iterator ipm, int const *nhp, int const *aszp);
	void _updatePrior(vector<bool> const &init, bool oneround = false){};

	double readInput(MYDATA &mydata, char const *input, vector<string> const &fS, vector<string> const &fA, double error, int cut, bool likelihood);
	void calPoislp(vector<SNPINFO> const &snpinfo, vector<vector<int> > const &counts, double error, double lambda, int cut, double Sdensity, double *poisp, double *nullp, double *errp, int aN);
	void outputHaplotype(char const *fname, MYDATA const &mydata);
};

