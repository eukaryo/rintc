/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa

The namespace parasor_param{} is modified from:

Kawaguchi R. et al. (2016) Parallel computation of genome-scale RNA secondary structure
to detect structural constraints on human genome. BMC Bioinformatics, 17:203.
*/

#include"parameter.h"


namespace rintdwr {

namespace parasor_param {

static constexpr double K0 = 273.15;
static constexpr double GASCONST = 1.98717;//Cal/(K*mol)
static constexpr double Tmeasure = K0 + 37.0;//Turner2004は37度で、Turner2004しか考えてないので決め打ち
static constexpr double SCALE = 10.0;
static constexpr int INTINF = 1000000;
static constexpr double INF = std::numeric_limits<double>::max() / 100.0;
static constexpr int PARAM_MAXLOOP = 30;

static const int8_t BP_pair[8][8] = {
	///*  _  A  C  G  U  X  K  I */
	{ 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 5, 0, 0, 0 },
	{ 0, 0, 0, 1, 0, 0, 0, 0 },
	{ 0, 0, 2, 0, 3, 0, 0, 0 },
	{ 0, 6, 0, 4, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0 },
	{ 0, 0, 0, 0, 0, 0, 0, 0 }
};

static const int8_t rtype[7] = { 0, 2, 1, 4, 3, 6, 5 };

static const int8_t typetable[256] = {
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,

	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,1/*'A'==65*/,0,2/*'C'==67*/,0,0,
	0,3/*'G'==71*/,0,0,0,0,0,0,0,0,
	0,0,0,0,0,4/*'U'==85*/,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,

	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,

	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,

	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,

	0,0,0,0,0,0
};

int GetBaseType(const char c) {
	return typetable[c];

	//switch (c) {
	//case 'A':return 1;
	//case 'C':return 2;
	//case 'G':return 3;
	//case 'U':return 4;
	//default:
	//	assert(0);
	//	break;
	//}
	//return 0;

}
int GetPairType(const char c1, const char c2) {
	return BP_pair[typetable[c1]][typetable[c2]];

	//if (c1 == 'C'&&c2 == 'G')return 1;
	//if (c1 == 'G'&&c2 == 'C')return 2;
	//if (c1 == 'G'&&c2 == 'U')return 3;
	//if (c1 == 'U'&&c2 == 'G')return 4;
	//if (c1 == 'A'&&c2 == 'U')return 5;
	//if (c1 == 'U'&&c2 == 'A')return 6;
	//return 0;

}
int GetPairTypeReverse(const char c1, const char c2) {
	return rtype[BP_pair[typetable[c1]][typetable[c2]]];

	//if (c1 == 'C'&&c2 == 'G')return 2;
	//if (c1 == 'G'&&c2 == 'C')return 1;
	//if (c1 == 'G'&&c2 == 'U')return 4;
	//if (c1 == 'U'&&c2 == 'G')return 3;
	//if (c1 == 'A'&&c2 == 'U')return 6;
	//if (c1 == 'U'&&c2 == 'A')return 5;
	//return 0;

}

static double kT;
static double TT;
static double loghairpin[PARAM_MAXLOOP + 1];
static double logmismatchH[7][5][5];
static double logmismatchI[7][5][5];
static double logmismatchM[7][5][5];
static double logmismatch1nI[7][5][5];
static double logmismatch23I[7][5][5];
static double logmismatchExt[7][5][5];
static double Triloop[40];
static double Tetraloop[40];
static double Hexaloop[40];
static double logstack[7][7];
static double logbulge[PARAM_MAXLOOP + 1];

static double logint11[8][8][5][5];
static double logint21[8][8][5][5][5];
static double logint22[8][8][5][5][5][5];
static double loginternal[PARAM_MAXLOOP + 1];
static double logdangle5[8][5];
static double logdangle3[8][5];
static double logninio[PARAM_MAXLOOP + 1];
static double logMLintern;
static double logMLclosing;
static double logML_BASE;

//ENTHALPIES
static double logstack_enthalpies[7][7];
static double logmismatchH_enthalpies[7][5][5];
static double logmismatchI_enthalpies[7][5][5];
static double logmismatch1nI_enthalpies[7][5][5];
static double logmismatch23I_enthalpies[7][5][5];
static double logmismatchM_enthalpies[7][5][5];
static double logmismatchExt_enthalpies[7][5][5];
static double logdangle5_enthalpies[8][5];
static double logdangle3_enthalpies[8][5];
static double logint11_enthalpies[8][8][5][5];
static double logint21_enthalpies[8][8][5][5][5];
static double logint22_enthalpies[8][8][5][5][5][5];
static double loghairpin_enthalpies[PARAM_MAXLOOP + 1];
static double logbulge_enthalpies[PARAM_MAXLOOP + 1];
static double loginternal_enthalpies[PARAM_MAXLOOP + 1];
static double logMLintern_enthalpy;
static double logMLclosing_enthalpy;
static double logML_BASE_enthalpy;
static double Triloop_enthalpies[40];
static double Tetraloop_enthalpies[40];
static double Hexaloop_enthalpies[40];

static double F_ninio37, F_ninio_enthalpy, MAX_NINIO;
static double TermAU_energy, TermAU_enthalpy;

//params, log scale
static double parhairpin[PARAM_MAXLOOP + 1];
static double parmismatchH[7][5][5];
static double parmismatchI[7][5][5];
static double parmismatchM[7][5][5];
static double parmismatch1nI[7][5][5];
static double parmismatch23I[7][5][5];
static double parmismatchExt[7][5][5];
static double parTriloop[40];
static double parTetraloop[40];
static double parHexaloop[40];
static double parstack[7][7];
static double parbulge[PARAM_MAXLOOP + 1];

static double parint11[8][8][5][5];
static double parint21[8][8][5][5][5];
static double parint22[8][8][5][5][5][5];
static double parinternal[PARAM_MAXLOOP + 1];
static double pardangle5[8][5];
static double pardangle3[8][5];
static double parninio[PARAM_MAXLOOP + 1];
static double parMLintern;
static double parMLclosing;
static double parML_BASE;
static double parTermAU;

static char Triloops[400];
static char Tetraloops[400];
static char Hexaloops[400];

static const bool no_closingGU = false;
static const bool tetra = true;

static double lxc37 = 107.856;
bool initialized = false;
bool inittermau = false;
bool old_param = false;

bool counting = false;

static double SMOOTH(const double X) {
	if ((X / SCALE) < -1.2283697)return 0.0;
	if ((X / SCALE) > 0.8660254)return X;
	return SCALE*0.38490018*(sin((X) / SCALE - 0.34242663) + 1.0)*(sin((X) / SCALE - 0.34242663) + 1.0);
}
static double SmoothingEnergy(const double enthalpy, const double energy37) {
	if (energy37 == INTINF) return 0.0;
	const double GT = enthalpy - (enthalpy - energy37) * TT;
	return SMOOTH(-GT) * 10.0 / kT;
}
static double Energy(const double enthalpy, const double energy37) {
	if (energy37 == INTINF) return 0.0;
	const double GT = enthalpy - (enthalpy - energy37) * TT;
	return -GT * 10.0 / kT;
}

static void ComputeParams() {
	for (int i = 0; i <= PARAM_MAXLOOP; ++i) {
		parhairpin[i] = Energy(loghairpin_enthalpies[i], loghairpin[i]);
		parbulge[i] = Energy(logbulge_enthalpies[i], logbulge[i]);
		parinternal[i] = Energy(loginternal_enthalpies[i], loginternal[i]);

		const double ninio_GT = (F_ninio_enthalpy - (F_ninio_enthalpy - F_ninio37) * TT);
		parninio[i] = -std::min(MAX_NINIO, i * ninio_GT) * 10.0 / kT;
	}
	for (int i = 1; i <= 6; ++i) {
		for (int j = 1; j <= 6; ++j) {
			parstack[i][j] = Energy(logstack_enthalpies[i][j], logstack[i][j]);
			for (int k = 1; k <= 4; ++k) {
				for (int l = 1; l <= 4; ++l) {
					parint11[i][j][k][l] = Energy(logint11_enthalpies[i][j][k][l], logint11[i][j][k][l]);
					for (int m = 1; m <= 4; ++m) {
						parint21[i][j][k][l][m] = Energy(logint21_enthalpies[i][j][k][l][m], logint21[i][j][k][l][m]);
						for (int n = 1; n <= 4; ++n) {
							parint22[i][j][k][l][m][n] = Energy(logint22_enthalpies[i][j][k][l][m][n], logint22[i][j][k][l][m][n]);
						}
					}
				}
			}
		}
		for (int j = 1; j <= 4; ++j) {
			pardangle5[i][j] = SmoothingEnergy(logdangle5_enthalpies[i][j], logdangle5[i][j]);
			pardangle3[i][j] = SmoothingEnergy(logdangle3_enthalpies[i][j], logdangle3[i][j]);
			for (int k = 1; k <= 4; ++k) {
				parmismatchH[i][j][k] = Energy(logmismatchH_enthalpies[i][j][k], logmismatchH[i][j][k]);
				parmismatchI[i][j][k] = Energy(logmismatchI_enthalpies[i][j][k], logmismatchI[i][j][k]);
				parmismatchM[i][j][k] = SmoothingEnergy(logmismatchM_enthalpies[i][j][k], logmismatchM[i][j][k]);
				parmismatch1nI[i][j][k] = Energy(logmismatch1nI_enthalpies[i][j][k], logmismatch1nI[i][j][k]);
				parmismatch23I[i][j][k] = Energy(logmismatch23I_enthalpies[i][j][k], logmismatch23I[i][j][k]);
				parmismatchExt[i][j][k] = SmoothingEnergy(logmismatchExt_enthalpies[i][j][k], logmismatchExt[i][j][k]);
			}
		}
	}
	parMLintern = Energy(logMLintern_enthalpy, logMLintern);
	parMLclosing = Energy(logMLclosing_enthalpy, logMLclosing);
	parML_BASE = Energy(logML_BASE_enthalpy, logML_BASE);
	parTermAU = Energy(TermAU_enthalpy, TermAU_energy);
	for (int i = 0; i < 40; ++i) {
		parTriloop[i] = Energy(Triloop_enthalpies[i], Triloop[i]);
		parTetraloop[i] = Energy(Tetraloop_enthalpies[i], Tetraloop[i]);
		//parTetraloop[i] = Energy(-400, Tetraloop[i]);
		parHexaloop[i] = Energy(Hexaloop_enthalpies[i], Hexaloop[i]);
	}
}

class Convert {
private:
	std::ifstream ifs;
	void GetWords(std::string& str, std::vector<std::string>& words) {
		words.clear();
		std::istringstream iss(str);
		copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(), back_inserter(words));
	}
	int param_type(std::string& str) {
		if (str.length() == 0 || str[0] != '#') return -1;
		std::vector<std::string> words;
		GetWords(str, words);
		if (words.size() <= 1) return -1;
		std::string id = words[1];
		if (id == "stack" || id == "stack_energies") return Stack;
		else if (id == "hairpin") return Hairpin;
		else if (id == "bulge") return Bulge;
		else if (id == "interior" || id == "internal_loop") return Interior;
		else if (id == "mismatch_exterior") return MisE;
		else if (id == "mismatch_hairpin") return MisH;
		else if (id == "mismatch_interior") return MisI;
		else if (id == "mismatch_interior_1n") return Mis1n;
		else if (id == "mismatch_interior_23") return MisI23;
		else if (id == "mismatch_multi") return MisM;
		else if (id == "int11" || id == "int11_energies") return Int11;
		else if (id == "int21" || id == "int21_energies") return Int21;
		else if (id == "int22" || id == "int22_energies") return Int22;
		else if (id == "dangle5") return Dan5;
		else if (id == "dangle3") return Dan3;
		else if (id == "ML_params") return ML;
		else if (id == "NINIO") return Ninio;
		else if (id == "Triloops") return Tri;
		else if (id == "Tetraloops") return Tetra;
		else if (id == "Hexaloops") return Hexa;
		else if (id == "Misc") return Misc;

		//ENTHALPIES
		else if (id == "stack_enthalpies") return Stack_ent;
		else if (id == "mismatch_hairpin_enthalpies") return MisH_ent;
		else if (id == "mismatch_interior_enthalpies") return MisI_ent;
		else if (id == "mismatch_interior_1n_enthalpies") return Mis1n_ent;
		else if (id == "mismatch_interior_23_enthalpies") return MisI23_ent;
		else if (id == "mismatch_multi_enthalpies") return MisM_ent;
		else if (id == "mismatch_exterior_enthalpies") return MisE_ent;
		else if (id == "dangle5_enthalpies") return Dan5_ent;
		else if (id == "dangle3_enthalpies") return Dan3_ent;
		else if (id == "int11_enthalpies" || id == "int11_energies_enthalpies") return Int11_ent;
		else if (id == "int21_enthalpies" || id == "int21_energies_enthalpies") return Int21_ent;
		else if (id == "int22_enthalpies" || id == "int22_energies_enthalpies") return Int22_ent;
		else if (id == "hairpin_enthalpies") return Hairpin_ent;
		else if (id == "bulge_enthalpies") return Bulge_ent;
		else if (id == "interior_enthalpies" || id == "internal_loop_enthalpies") return Interior_ent;

		return -1;
	}
	void GetArray(double *arr_, const int size, const bool smooth = false) {
		std::vector<std::string> words;
		for (int i = 0; i < size; ) {
			std::string str;
			if (!getline(ifs, str)) break;
			if (str.size() > 0 && str[str.size() - 1] == '\r') {
				str.erase(str.size() - 1);
			}
			else if (str.length() < 2) break;
			GetWords(str, words);
			int prev = i;
			for (; i < size && i - prev < (int)words.size(); i++) {
				if (words[i - prev].find("/*") != std::string::npos) break;
				else {
					int temp;
					if (words[i - prev] == "INF") temp = INTINF;
					else if (words[i - prev] == "DEF") temp = -50;
					else temp = atoi(words[i - prev].c_str());
					arr_[i] = double(temp);
				}
			}
		}
	}
	void Read1Dim(double *arr_, const int dim, const int shift, const int post = 0) {
		GetArray(arr_ + shift, dim - shift - post);
	}
	void Read2Dim(double* arr_, const int dim1, const int dim2, const int shift1, const int shift2,
		const int post1 = 0, const int post2 = 0)
	{
		if (shift1 + shift2 == 0 && post1 + post2 == 0)
			Read1Dim(arr_, dim1*dim2, 0);
		else {
			for (int i = shift1; i < dim1 - post1; i++)
				Read1Dim(arr_ + (i*dim2), dim2, shift2, post2);
		}
	}
	void Read3Dim(double *arr_, const int dim1, const int dim2, const int dim3, const int shift1, const int shift2, const int shift3,
		const int post1 = 0, const int post2 = 0, const int post3 = 0)
	{
		if (shift1 + shift2 + shift3 == 0 && post1 + post2 + post3 == 0)
			Read1Dim(arr_, dim1*dim2*dim3, 0);
		else {
			for (int i = shift1; i < dim1 - post1; i++) {
				Read2Dim(arr_ + (i*dim2*dim3), dim2, dim3, shift2, shift3, post2, post3);
			}
		}
	}
	void Read4Dim(double *arr_, const int dim1, const int dim2, const int dim3, const int dim4,
		const int shift1, const int shift2, const int shift3, const int shift4,
		const int post1 = 0, const int post2 = 0, const int post3 = 0, const int post4 = 0)
	{
		if (shift1 + shift2 + shift3 + shift4 == 0 && post1 + post2 + post3 + post4 == 0)
			Read1Dim(arr_, dim1*dim2*dim3*dim4, 0);
		else {
			for (int i = shift1; i < dim1 - post1; i++) {
				Read3Dim(arr_ + (i*dim2*dim3*dim4), dim2, dim3, dim4, shift2, shift3, shift4,
					post2, post3, post4);
			}
		}
	}
	void Read5Dim(double *arr_, const int dim1, const int dim2, const int dim3, const int dim4, const int dim5,
		const int shift1, const int shift2, const int shift3, const int shift4, const int shift5,
		const int post1 = 0, const int post2 = 0, const int post3 = 0, const int post4 = 0, const int post5 = 0)
	{
		if (shift1 + shift2 + shift3 + shift4 + shift5 == 0 && post1 + post2 + post3 + post4 + post5 == 0)
			Read1Dim(arr_, dim1*dim2*dim3*dim4*dim5, 0);
		else {
			for (int i = shift1; i < dim1 - post1; i++) {
				Read4Dim(arr_ + (i*dim2*dim3*dim4*dim5), dim2, dim3, dim4, dim5,
					shift2, shift3, shift4, shift5,
					post2, post3, post4, post5);
			}
		}
	}
	void Read6Dim(double *arr_, const int dim1, const int dim2, const int dim3, const int dim4, const int dim5, const int dim6,
		const int shift1, const int shift2, const int shift3, const int shift4, const int shift5, const int shift6,
		const int post1 = 0, const int post2 = 0, const int post3 = 0, const int post4 = 0, const int post5 = 0, const int post6 = 0)
	{
		if (shift1 + shift2 + shift3 + shift4 + shift5 + shift6 == 0 && post1 + post2 + post3 + post4 + post5 + post6 == 0)
			Read1Dim(arr_, dim1*dim2*dim3*dim4*dim5*dim6, 0);
		else {
			for (int i = shift1; i < dim1 - post1; i++) {
				Read5Dim(arr_ + (i*dim2*dim3*dim4*dim5*dim6), dim2, dim3, dim4, dim5, dim6,
					shift2, shift3, shift4, shift5, shift6,
					post2, post3, post4, post5, post6);
			}
		}
	}
	void Readninio()
	{
		std::string str;
		std::vector<std::string> words, pwords;
		//int F_ninio37, MAX_NINIO;
		while (getline(ifs, str)) {
			if (str.size() > 0 && str[str.size() - 1] == '\r') {
				str.erase(str.size() - 1);
			}
			if (str == "") break;
			GetWords(str, words);
			if (str.find("*") != std::string::npos) {
				pwords = words;
				continue;
			}
			if ((int)pwords.size() == (int)words.size() + 2) {

				//cout << "Turner2004 only" << endl;
				//assert(0);

				assert((int)words.size() >= 2);
				for (int i = 0; i < (int)words.size(); i++) {
					int value = atoi(words[i].c_str());
					if (pwords[i + 1] == "m") F_ninio37 = value;
					else if (pwords[i + 1] == "m_dH") F_ninio_enthalpy = value;
					else if (pwords[i + 1] == "max") MAX_NINIO = value;
				}
			}
			else {
				assert((int)words.size() > 2);
				F_ninio37 = atoi(words[0].c_str());
				MAX_NINIO = atoi(words[2].c_str());

				F_ninio_enthalpy = atoi(words[1].c_str());
			}
			break;
		}
	}
	void ReadML()
	{
		std::string str;
		std::vector<std::string> words, pwords;
		while (getline(ifs, str)) {
			if (str.size() > 0 && str[str.size() - 1] == '\r') {
				str.erase(str.size() - 1);
			}
			if (str == "") break;
			GetWords(str, words);
			if (str.find("*") != std::string::npos) {
				pwords = words;
				continue;
			}

			if ((int)pwords.size() == (int)words.size() + 2) {

				//cout << "Turner2004 only" << endl;
				//assert(0);

				assert((int)words.size() > 3);
				for (int i = 0; i < (int)words.size(); i++) {
					double value = atof(words[i].c_str());
					if (pwords[i + 1] == "cu") logML_BASE = value;
					else if (pwords[i + 1] == "cc") logMLclosing = value;
					else if (pwords[i + 1] == "ci") logMLintern = value;
					else if (pwords[i + 1] == "cu_dH") logML_BASE_enthalpy = value;
					else if (pwords[i + 1] == "cc_dH") logMLclosing_enthalpy = value;
					else if (pwords[i + 1] == "ci_dH") logMLintern_enthalpy = value;
					else if (pwords[i + 1] == "TerminalAU") SetAU(value);
					// else cout << "Unidentified parameter: " << words[i] << " " << pwords[i+1] << endl;
				}
			}
			else {
				assert((int)words.size() > 5);
				logML_BASE = atof(words[0].c_str());
				logMLclosing = atof(words[2].c_str());
				logMLintern = atof(words[4].c_str());

				logML_BASE_enthalpy = atof(words[1].c_str());
				logMLclosing_enthalpy = atof(words[3].c_str());
				logMLintern_enthalpy = atof(words[5].c_str());
			}
			break;
		}
	}
	void ReadMisc(const bool lxc = false)
	{
		std::string str;
		std::vector<std::string> words, pwords;
		while (getline(ifs, str)) {
			if (str.size() > 0 && str[str.size() - 1] == '\r') {
				str.erase(str.size() - 1);
			}
			if (str == "") break;
			GetWords(str, words);
			if (str.find("*") != std::string::npos) continue;
			if ((int)pwords.size() == (int)words.size() + 2) {

				std::cout << "Turner2004 only" << std::endl;
				assert(0);

				assert((int)words.size() > 2);
				for (int i = 0; i < (int)words.size(); i++) {
					if (lxc) {
						if (pwords[i + 1] == "LXC" || pwords[i + 1] == "lxc") lxc37 = atof(words[i].c_str());
					}
					else if (pwords[i + 1] == "TerminalAU") SetAU(words[i]);
				}
			}
			else {
				assert((int)words.size() >= 4);
				if (lxc) {
					if ((int)words.size() > 4) {
						lxc37 = atof(words[4].c_str());
					}
				}
				//else SetAU(words[2]);

				TermAU_energy = atof(words[2].c_str());
				TermAU_enthalpy = atof(words[3].c_str());

			}
			break;
		}
	}
	void ReadString(double* arr_, double* arr_enthalpies, std::string & loopstr)
	{
		std::string str;
		std::vector<std::string> words;
		for (int i = 0; getline(ifs, str); i++) {
			if (str.size() > 0 && str[str.size() - 1] == '\r') {
				str.erase(str.size() - 1);
			}
			if (str == "") break;
			else if (str.find("*") != std::string::npos) {
				i--; continue;
			}
			GetWords(str, words);
			assert((int)words.size() > 1);
			loopstr += words[0] + " ";
			arr_[i] = atof(words[1].c_str());
			arr_enthalpies[i] = atof(words[2].c_str());
		}
	}
	void Read2DimSmooth(double *arr_, const int dim1, const int dim2, const int shift1, const int shift2, const int post1 = 0, const int post2 = 0)
	{
		if (shift1 + shift2 == 0 && post1 + post2 == 0)
			GetArray(arr_, dim1*dim2, true);
		else {
			for (int i = shift1; i < dim1 - post1; i++)
				GetArray(arr_ + (i*dim2) + shift2, dim2 - shift2 - post2, true);
		}
	}
	void Read3DimSmooth(double *arr_, const int dim1, const int dim2, const int dim3, const int shift1, const int shift2, const int shift3,
		const int post1 = 0, const int post2 = 0, const int post3 = 0)
	{
		if (shift1 + shift2 + shift3 == 0 && post1 + post2 + post3 == 0)
			GetArray(arr_, dim1*dim2*dim3, true);
		else {
			for (int i = shift1; i < dim1 - post1; i++) {
				Read2DimSmooth(arr_ + (i*dim2*dim3), dim2, dim3, shift2, shift3, post2, post3);
			}
		}
	}
	void FillINF()
	{
		std::fill(&(logstack[0][0]), &(logstack[0][0]) + 7 * 7, -INF);
		std::fill(&(logmismatchH[0][0][0]), &(logmismatchH[0][0][0]) + 7 * 5 * 5, -INF);
		std::fill(&(logmismatchI[0][0][0]), &(logmismatchI[0][0][0]) + 7 * 5 * 5, -INF);
		std::fill(&(logmismatch1nI[0][0][0]), &(logmismatch1nI[0][0][0]) + 7 * 5 * 5, -INF);
		std::fill(&(logmismatch23I[0][0][0]), &(logmismatch23I[0][0][0]) + 7 * 5 * 5, -INF);
		std::fill(&(logmismatchM[0][0][0]), &(logmismatchM[0][0][0]) + 7 * 5 * 5, -INF);
		std::fill(&(logmismatchExt[0][0][0]), &(logmismatchExt[0][0][0]) + 7 * 5 * 5, -INF);
		std::fill(&(logdangle5[0][0]), &(logdangle5[0][0]) + 8 * 5, -INF);
		std::fill(&(logdangle3[0][0]), &(logdangle3[0][0]) + 8 * 5, -INF);
		std::fill(&(logint11[0][0][0][0]), &(logint11[0][0][0][0]) + 8 * 8 * 5 * 5, -INF);
		std::fill(&(logint21[0][0][0][0][0]), &(logint21[0][0][0][0][0]) + 8 * 8 * 5 * 5 * 5, -INF);
		std::fill(&(logint22[0][0][0][0][0][0]), &(logint22[0][0][0][0][0][0]) + 8 * 8 * 5 * 5 * 5 * 5, -INF);//ParasoRのバグ
		std::fill(&(loghairpin[0]), &(loghairpin[0]) + PARAM_MAXLOOP + 1, -INF);
		std::fill(&(logbulge[0]), &(logbulge[0]) + PARAM_MAXLOOP + 1, -INF);
		std::fill(&(loginternal[0]), &(loginternal[0]) + PARAM_MAXLOOP + 1, -INF);
		std::fill(&(logninio)[0], &(logninio[0]) + PARAM_MAXLOOP + 1, -INF);
		std::fill(&(Triloop[0]), &(Triloop[0]) + 40, -INF);
		std::fill(&(Tetraloop[0]), &(Tetraloop[0]) + 40, -INF);
		std::fill(&(Hexaloop[0]), &(Hexaloop[0]) + 40, -INF);

		//ENTHALPIES
		std::fill(&(logstack_enthalpies[0][0]), &(logstack_enthalpies[0][0]) + 7 * 7, -INF);
		std::fill(&(logmismatchH_enthalpies[0][0][0]), &(logmismatchH_enthalpies[0][0][0]) + 7 * 5 * 5, -INF);
		std::fill(&(logmismatchI_enthalpies[0][0][0]), &(logmismatchI_enthalpies[0][0][0]) + 7 * 5 * 5, -INF);
		std::fill(&(logmismatch1nI_enthalpies[0][0][0]), &(logmismatch1nI_enthalpies[0][0][0]) + 7 * 5 * 5, -INF);
		std::fill(&(logmismatch23I_enthalpies[0][0][0]), &(logmismatch23I_enthalpies[0][0][0]) + 7 * 5 * 5, -INF);
		std::fill(&(logmismatchM_enthalpies[0][0][0]), &(logmismatchM_enthalpies[0][0][0]) + 7 * 5 * 5, -INF);
		std::fill(&(logmismatchExt_enthalpies[0][0][0]), &(logmismatchExt_enthalpies[0][0][0]) + 7 * 5 * 5, -INF);
		std::fill(&(logdangle5_enthalpies[0][0]), &(logdangle5_enthalpies[0][0]) + 8 * 5, -INF);
		std::fill(&(logdangle3_enthalpies[0][0]), &(logdangle3_enthalpies[0][0]) + 8 * 5, -INF);
		std::fill(&(logint11_enthalpies[0][0][0][0]), &(logint11_enthalpies[0][0][0][0]) + 8 * 8 * 5 * 5, -INF);
		std::fill(&(logint21_enthalpies[0][0][0][0][0]), &(logint21_enthalpies[0][0][0][0][0]) + 8 * 8 * 5 * 5 * 5, -INF);
		std::fill(&(logint22_enthalpies[0][0][0][0][0][0]), &(logint22_enthalpies[0][0][0][0][0][0]) + 8 * 8 * 5 * 5 * 5 * 5, -INF);
		std::fill(&(loghairpin_enthalpies[0]), &(loghairpin_enthalpies[0]) + PARAM_MAXLOOP + 1, -INF);
		std::fill(&(logbulge_enthalpies[0]), &(logbulge_enthalpies[0]) + PARAM_MAXLOOP + 1, -INF);
		std::fill(&(loginternal_enthalpies[0]), &(loginternal_enthalpies[0]) + PARAM_MAXLOOP + 1, -INF);
		std::fill(&(Triloop_enthalpies[0]), &(Triloop_enthalpies[0]) + 40, -INF);
		std::fill(&(Tetraloop_enthalpies[0]), &(Tetraloop_enthalpies[0]) + 40, -INF);
		std::fill(&(Hexaloop_enthalpies[0]), &(Hexaloop_enthalpies[0]) + 40, -INF);
	}
	void ReadOnlyMisc(const std::string& file)
	{
		std::string str;
		ifs.open(file.c_str());
		while (getline(ifs, str)) {
			if (str.size() > 0 && str[str.size() - 1] == '\r') {
				str.erase(str.size() - 1);
			}
			int num = param_type(str);
			if (num == Misc) {
				ReadMisc(true);
				break;
			}
		}
		ifs.close();
		FillINF();
	}
	void SetAU(const std::string value)
	{
		assert(0);
		//logTermAU = LogEnergy(atoi(value.c_str()));
		//inittermau = true;
	}
	void SetAU(const double value)
	{
		assert(0);
		//logTermAU = value;
		//inittermau = true;
	}
public:
	enum Params {
		Stack, MisH, MisI, Mis1n, MisI23, MisM, MisE, Dan5, Dan3, Int11, Int21, Int22,
		Hairpin, Bulge, Interior, Ninio, ML, Misc, Tri, Tetra, Hexa,
		Stack_ent, MisH_ent, MisI_ent, Mis1n_ent, MisI23_ent, MisM_ent, MisE_ent, Dan5_ent, Dan3_ent, Int11_ent, Int21_ent, Int22_ent,
		Hairpin_ent, Bulge_ent, Interior_ent
	};//ENTHALPIES
	Convert() {}
	virtual ~Convert() {}
	bool ConvertParamFile(const std::string& file)
	{
		std::string str;
		ReadOnlyMisc(file);
		ifs.open(file.c_str());
		old_param = true;
		if (getline(ifs, str) && str.find("## RNAfold parameter file v2.0") != std::string::npos)
			old_param = false;
		while (getline(ifs, str)) {
			if (str.size() > 0 && str[str.size() - 1] == '\r') {
				str.erase(str.size() - 1);
			}
			int num = param_type(str);
			if (num == Stack) {
				if (old_param) Read2Dim(&(logstack[0][0]), 7, 7, 0, 0);
				else Read2Dim(&(logstack[0][0]), 7, 7, 1, 1);
			}
			else if (num == MisH) {
				if (old_param) Read3Dim(&(logmismatchH[0][0][0]), 7, 5, 5, 0, 0, 0);
				else Read3Dim(&(logmismatchH[0][0][0]), 7, 5, 5, 1, 0, 0);
			}
			else if (num == MisI) {
				if (old_param) Read3Dim(&(logmismatchI[0][0][0]), 7, 5, 5, 0, 0, 0);
				else Read3Dim(&(logmismatchI[0][0][0]), 7, 5, 5, 1, 0, 0);
			}
			else if (num == Mis1n) {
				Read3Dim(&(logmismatch1nI[0][0][0]), 7, 5, 5, 1, 0, 0);
			}
			else if (num == MisI23) {
				Read3Dim(&(logmismatch23I[0][0][0]), 7, 5, 5, 1, 0, 0);
			}
			else if (num == MisM) {
				Read3DimSmooth(&(logmismatchM[0][0][0]), 7, 5, 5, 1, 0, 0);
			}
			else if (num == MisE) {
				Read3DimSmooth(&(logmismatchExt[0][0][0]), 7, 5, 5, 1, 0, 0);
			}
			else if (num == Dan5) {
				if (old_param) Read2DimSmooth(&(logdangle5[0][0]), 8, 5, 0, 0);
				else Read2DimSmooth(&(logdangle5[0][0]), 8, 5, 1, 0);
			}
			else if (num == Dan3) {
				if (old_param) Read2DimSmooth(&(logdangle3[0][0]), 8, 5, 0, 0);
				else Read2DimSmooth(&(logdangle3[0][0]), 8, 5, 1, 0);
			}
			else if (num == Int11) {
				Read4Dim(&(logint11[0][0][0][0]), 8, 8, 5, 5,
					1, 1, 0, 0);
			}
			else if (num == Int21) {
				Read5Dim(&(logint21[0][0][0][0][0]), 8, 8, 5, 5, 5,
					1, 1, 0, 0, 0);
			}
			else if (num == Int22) {
				if (old_param) {
					Read6Dim(&(logint22[0][0][0][0][0][0]), 8, 8, 5, 5, 5, 5,
						1, 1, 1, 1, 1, 1,
						0, 0, 0, 0, 0, 0);
				}
				else {
					Read6Dim(&(logint22[0][0][0][0][0][0]), 8, 8, 5, 5, 5, 5,
						1, 1, 1, 1, 1, 1,
						1, 1, 0, 0, 0, 0);
				}
			}
			else if (num == Hairpin) {
				Read1Dim(&(loghairpin[0]), PARAM_MAXLOOP + 1, 0);
			}
			else if (num == Bulge) {
				Read1Dim(&(logbulge[0]), PARAM_MAXLOOP + 1, 0);
			}
			else if (num == Interior) {
				Read1Dim(&(loginternal[0]), PARAM_MAXLOOP + 1, 0);
			}
			else if (num == Ninio) {
				Readninio();
			}
			else if (num == ML) {
				ReadML();
			}
			else if (num == Misc) {
				ReadMisc();
			}
			else if (num == Tri) {
				std::string tmp;
				ReadString(&(Triloop[0]), &(Triloop_enthalpies[0]), tmp);
				try {
					//strcpy_s(Triloops, tmp.c_str());
					for (int a = 0; a < 400; ++a)Triloops[a] = '\0';
					for (int a = 0; a < tmp.size(); ++a)Triloops[a] = tmp[a];
				}
				catch (...) {
					assert(0);
				}
			}
			else if (num == Tetra) {
				std::string tmp;
				ReadString(&(Tetraloop[0]), &(Tetraloop_enthalpies[0]), tmp);
				try {
					//strcpy_s(Tetraloops, tmp.c_str());
					for (int a = 0; a < 400; ++a)Tetraloops[a] = '\0';
					for (int a = 0; a < tmp.size(); ++a)Tetraloops[a] = tmp[a];
				}
				catch (...) {
					assert(0);
				}
			}
			else if (num == Hexa) {
				std::string tmp;
				ReadString(&(Hexaloop[0]), &(Hexaloop_enthalpies[0]), tmp);
				try {
					//strcpy_s(Hexaloops, tmp.c_str());
					for (int a = 0; a < 400; ++a)Hexaloops[a] = '\0';
					for (int a = 0; a < tmp.size(); ++a)Hexaloops[a] = tmp[a];
				}
				catch (...) {
					assert(0);
				}
			}

			//ENTHALPIES
			else if (num == Stack_ent) {
				if (old_param) Read2Dim(&(logstack_enthalpies[0][0]), 7, 7, 0, 0);
				else Read2Dim(&(logstack_enthalpies[0][0]), 7, 7, 1, 1);
			}
			else if (num == MisH_ent) {
				if (old_param) Read3Dim(&(logmismatchH_enthalpies[0][0][0]), 7, 5, 5, 0, 0, 0);
				else Read3Dim(&(logmismatchH_enthalpies[0][0][0]), 7, 5, 5, 1, 0, 0);
			}
			else if (num == MisI_ent) {
				if (old_param) Read3Dim(&(logmismatchI_enthalpies[0][0][0]), 7, 5, 5, 0, 0, 0);
				else Read3Dim(&(logmismatchI_enthalpies[0][0][0]), 7, 5, 5, 1, 0, 0);
			}
			else if (num == Mis1n_ent) {
				Read3Dim(&(logmismatch1nI_enthalpies[0][0][0]), 7, 5, 5, 1, 0, 0);
			}
			else if (num == MisI23_ent) {
				Read3Dim(&(logmismatch23I_enthalpies[0][0][0]), 7, 5, 5, 1, 0, 0);
			}
			else if (num == MisM_ent) {
				Read3DimSmooth(&(logmismatchM_enthalpies[0][0][0]), 7, 5, 5, 1, 0, 0);
			}
			else if (num == MisE_ent) {
				Read3DimSmooth(&(logmismatchExt_enthalpies[0][0][0]), 7, 5, 5, 1, 0, 0);
			}
			else if (num == Dan5_ent) {
				if (old_param) Read2DimSmooth(&(logdangle5_enthalpies[0][0]), 8, 5, 0, 0);
				else Read2DimSmooth(&(logdangle5_enthalpies[0][0]), 8, 5, 1, 0);
			}
			else if (num == Dan3_ent) {
				if (old_param) Read2DimSmooth(&(logdangle3_enthalpies[0][0]), 8, 5, 0, 0);
				else Read2DimSmooth(&(logdangle3_enthalpies[0][0]), 8, 5, 1, 0);
			}
			else if (num == Int11_ent) {
				Read4Dim(&(logint11_enthalpies[0][0][0][0]), 8, 8, 5, 5,
					1, 1, 0, 0);
			}
			else if (num == Int21_ent) {
				Read5Dim(&(logint21_enthalpies[0][0][0][0][0]), 8, 8, 5, 5, 5,
					1, 1, 0, 0, 0);
			}
			else if (num == Int22_ent) {
				if (old_param) {
					Read6Dim(&(logint22_enthalpies[0][0][0][0][0][0]), 8, 8, 5, 5, 5, 5,
						1, 1, 1, 1, 1, 1,
						0, 0, 0, 0, 0, 0);
				}
				else {
					Read6Dim(&(logint22_enthalpies[0][0][0][0][0][0]), 8, 8, 5, 5, 5, 5,
						1, 1, 1, 1, 1, 1,
						1, 1, 0, 0, 0, 0);
				}
			}
			else if (num == Hairpin_ent) {
				Read1Dim(&(loghairpin_enthalpies[0]), PARAM_MAXLOOP + 1, 0);
			}
			else if (num == Bulge_ent) {
				Read1Dim(&(logbulge_enthalpies[0]), PARAM_MAXLOOP + 1, 0);
			}
			else if (num == Interior_ent) {
				Read1Dim(&(loginternal_enthalpies[0]), PARAM_MAXLOOP + 1, 0);
			}


			else {
				int aaa = 0;
			}
		}
		ifs.close();
		return true;
	}
};

static inline bool IsAU(const int type) {
	return (type > 2);
}
static inline bool IsCloseGU(const int type) {
	return (type == 3 || type == 4);
}
static void ChangeEnergyParam(const std::string name = "")
{
	//bool doubt = false;
	std::string file;
	if (name.length() == 0 || name == "Turner2004"/* || name[0] == 'T'*/)
		//file = "energy_param/rna_turner2004_new.par";
		file = "energy_param/rna_turner2004.par";
	else if (name == "Andronescu" || name[0] == 'A')
		file = "energy_param/rna_andronescu2007.par";
	else if (name == "Turner1999")
		file = "energy_param/rna_turner1999.par";
	else {
		//file = name;
		//doubt = true;
		assert(0);
	}
	//cout << "#-Read Energy " << name << " -> " << file << endl;
	class Convert convert;
	if (!convert.ConvertParamFile(file)) {
		//cerr << "Format error: Energy param file " << file << endl;
		assert(0);
	}
	else {
		ComputeParams();
		initialized = true;
	}
	//if (doubt) PrintSummary();
}
static void SetTemperature(const double temp) {
	kT = (temp + K0) * GASCONST;
	TT = (temp + K0) / Tmeasure;
}

double ParMultiloopClosing() {
	if (counting)return 0.0;
	return parMLclosing;
}
double ParMultiloopInternal() {
	if (counting)return 0.0;
	return parMLintern;
}
double ParDangling(const int type, const int five, const int three, const bool ext, const std::string& sequence)
{
	if (counting)return 0.0;
	//five and three are 0-origin.
	//sequence = [ACGU]*

	assert(-1 <= five&&three <= sequence.size());

	double temp = 0.0;
	if (five >= 0 && three <= sequence.size() - 1 && !old_param) {
		if (ext) temp += parmismatchExt[type][GetBaseType(sequence[five])][GetBaseType(sequence[three])];
		else temp += parmismatchM[type][GetBaseType(sequence[five])][GetBaseType(sequence[three])];
		if (IsAU(type)) temp += parTermAU;
	}
	else {
		if (five >= 0) temp += pardangle5[type][GetBaseType(sequence[five])];
		if (three <= sequence.size() - 1) temp += pardangle3[type][GetBaseType(sequence[three])];
		if (IsAU(type)) temp += parTermAU;
	}
	// cout << "ext" << ext << " " << type << " " << five << "-" << three << " " <<sequence.seqget(five) << " " << sequence.seqget(three) << " " << temp << endl;
	return temp;
}
double ParHairpinEnergy(const int i, const int j, const std::string& sequence)
{
	if (counting)return 0.0;
	//i and j are 0-origin.
	//sequence = [ACGU]*

	const int type = GetPairType(sequence[i], sequence[j]);
	const int d = j - i - 1;
	double q = (d <= PARAM_MAXLOOP) ? parhairpin[d] : parhairpin[PARAM_MAXLOOP] - (lxc37*log(double(d) / double(PARAM_MAXLOOP))*10.0 / kT);
	if (d < 3) return q;
	if (tetra && d == 4) {
		const std::string sub_seq = sequence.substr(i, d + 2);
		const size_t tel = std::string(Tetraloops).find(sub_seq);
		if (tel != std::string::npos) {
			if (type != 7) return parTetraloop[tel / 7];
			else q += parTetraloop[tel / 7];
		}
	}
	if (tetra && d == 6) {
		const std::string sub_seq = sequence.substr(i, d + 2);
		const size_t tel = std::string(Hexaloops).find(sub_seq);
		if (tel != std::string::npos) return parHexaloop[tel / 9];
	}
	if (d == 3) {
		const std::string sub_seq = sequence.substr(i, d + 2);
		const size_t tel = std::string(Triloops).find(sub_seq);
		if (tel != std::string::npos) return parTriloop[tel / 6];
		if (IsAU(type)) q += parTermAU;
	}
	else {
		q += parmismatchH[type][GetBaseType(sequence[i + 1])][GetBaseType(sequence[j - 1])];
	}
	return q;
}
double ParLoopEnergy(const int i, const int j, const int p, const int q, const std::string& sequence)
{
	if (counting)return 0.0;
	assert(i < p&&p < q&&q < j);
	const int type1 = GetPairType(sequence[i], sequence[j]);
	const int type2 = GetPairTypeReverse(sequence[p], sequence[q]);
	const int u1 = p - i - 1;
	const int u2 = j - q - 1;
	const int u = std::max(u1, u2);
	assert(u1 + u2 <= PARAM_MAXLOOP);
	double z;
	if (u1 == 0 && u2 == 0) {
		return parstack[type1][type2];
	}
	else if (no_closingGU && (IsCloseGU(type1) || IsCloseGU(type2))) {
		return 0.0;
	}
	else if ((u1 == 0) || (u2 == 0)) { /* bulge */
		z = parbulge[u];
		if (u == 1) z += parstack[type1][type2];
		else {
			if (IsAU(type1)) z += parTermAU;
			if (IsAU(type2)) z += parTermAU;
		}
		return z;
	}
	else {
		if (u <= 2) {             /* short internal */
			if (u1 + u2 == 2) {
				z = parint11[type1][type2][GetBaseType(sequence[i + 1])][GetBaseType(sequence[j - 1])];
			}
			else if (u1 == 1 && u2 == 2) {
				z = parint21[type1][type2][GetBaseType(sequence[i + 1])][GetBaseType(sequence[q + 1])][GetBaseType(sequence[j - 1])];
			}
			else if (u1 == 2 && u2 == 1) {
				z = parint21[type2][type1][GetBaseType(sequence[q + 1])][GetBaseType(sequence[i + 1])][GetBaseType(sequence[p - 1])];
			}
			else {
				z = parint22[type1][type2][GetBaseType(sequence[i + 1])][GetBaseType(sequence[p - 1])][GetBaseType(sequence[q + 1])][GetBaseType(sequence[j - 1])];
			}
		}
		else {                  /* long internal */
			z = parinternal[u1 + u2];
			double temp1, temp2;
			if (u1 == 1 || u2 == 1) {
				temp1 = parmismatch1nI[type1][GetBaseType(sequence[i + 1])][GetBaseType(sequence[j - 1])];
				temp2 = parmismatch1nI[type2][GetBaseType(sequence[q + 1])][GetBaseType(sequence[p - 1])];
			}
			else if (u1 + u2 == 5) {
				temp1 = parmismatch23I[type1][GetBaseType(sequence[i + 1])][GetBaseType(sequence[j - 1])];
				temp2 = parmismatch23I[type2][GetBaseType(sequence[q + 1])][GetBaseType(sequence[p - 1])];
			}
			else {
				temp1 = parmismatchI[type1][GetBaseType(sequence[i + 1])][GetBaseType(sequence[j - 1])];
				temp2 = parmismatchI[type2][GetBaseType(sequence[q + 1])][GetBaseType(sequence[p - 1])];
			}
			const double temp3 = parninio[abs(u1 - u2)];
			z += temp1 + temp2 + temp3;
		}
	}
	return z;
}

void InitializeParameter(const std::string& parameter_file_name, const double temperature) {
	SetTemperature(temperature);
	ChangeEnergyParam(parameter_file_name);
	return;
}

}

}
