/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <stack>
#include <vector>
#include <map>
#include <set>
#include <cassert>
#include <iterator>
#include <functional>
#include <complex>
#include <random>
#include <chrono>
#include <sstream>
#include <thread>

#include <omp.h>

#include"complex_number.h"
#include"fourier_transform.h"
#include"misc.h"
#include"test.h"
#include"rintx.h"
#include"mccaskill_1990.h"
#include"sample_mccaskill.h"
#include"experiment.h"
#include"parameter.h"
#include"centroid_fold.h"

#ifdef _WIN64 
    #include<windows.h>
#endif

namespace rintdwr {

void TestAll() {

	const int num = 100;

	TestMaxHamming(num);
	TestSampling(num);
	TestSimpleMcCaskillWide(num);
	TestRintD1Dim(num);
	TestRintD2Dim(num);
	TestRintW1Dim(num);
	TestRintW2Dim(num);
	TestCentroidFold(num);
	TestHagioNonFourier(num);
	TestMaxHammingPK(num);
	TestRintD1DimPK(num);
}

int main_(int argc, char *argv[]) {

//#ifdef _WIN64 
//	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
//#endif


	//TestAll();

	if (argc == 2) {
		if (std::string(argv[1]) == std::string("test")) {
			TestAll();
			return 0;
		}
	}
	if (argc == 5) {

		const std::string sequence = std::string(argv[1]);
		const std::string structure = std::string(argv[2]);
		const int W = std::stoi(std::string(argv[3]));
		const std::string algo = std::string(argv[4]);

		RintX1DOptions options;
		options.temperature = 37.0;
		options.param_file_name = std::string("Turner2004");
		const int n = sequence.length();
		options.max_span = W;
		options.max_loop = n < 30 ? n : 30;
		options.sequence = sequence;
		options.reference_structure1 = structure;
		options.S1 = VerifyAndParseStructure(options.reference_structure1, options.sequence, options.max_span, options.max_loop);
		options.max_dim1 = ComputeMaxHammingDistance(options.sequence, options.S1, options.max_span, options.max_loop);

		if (algo == std::string("RintCwithDFT")) {
			typedef WideComplexNumber<Floating> Comp;
			options.allow_fft = false;
			const auto ans1 = ComputeRintW1Dim<Comp>(options);
			const auto ans2 = RintW1DimToBppm(ans1.first, ans1.second);
			for (int x = 0; x < ans2.first.size(); ++x) {
				std::cout << x << " " << ans2.first[x] << std::endl;
				if (ans2.first[x] == 0.0)continue;
				for (int i = 1; i <= n; ++i)for (int j = i + 1; j <= n && j - i <= options.max_span; ++j) {
					std::cout << x << " " << i << " " << j << " " << ans2.second[x][i][j - i] << std::endl;
				}
			}
			return 0;
		}
		if (algo == std::string("RintCwithFFT")) {
			typedef WideComplexNumber<Floating> Comp;
			options.allow_fft = true;
			const auto ans1 = ComputeRintW1Dim<Comp>(options);
			const auto ans2 = RintW1DimToBppm(ans1.first, ans1.second);
			for (int x = 0; x < ans2.first.size(); ++x) {
				std::cout << x << " " << ans2.first[x] << std::endl;
				if (ans2.first[x] == 0.0)continue;
				for (int i = 1; i <= n; ++i)for (int j = i + 1; j <= n && j - i <= options.max_span; ++j) {
					std::cout << x << " " << i << " " << j << " " << ans2.second[x][i][j - i] << std::endl;
				}
			}
			return 0;
		}
		else if (algo == std::string("RintCwithoutFourier")) {
			typedef WideComplexNumber<Floating> Comp;
			options.allow_fft = false;
			const auto ans1 = ComputeRintW1DimNonFourier<Comp>(options);
			const auto ans2 = RintW1DimToBppm(ans1.first, ans1.second);
			for (int x = 0; x < ans2.first.size(); ++x) {
				std::cout << x << " " << ans2.first[x] << std::endl;
				if (ans2.first[x] == 0.0)continue;
				for (int i = 1; i <= n; ++i)for (int j = i + 1; j <= n && j - i <= options.max_span; ++j) {
					std::cout << x << " " << i << " " << j << " " << ans2.second[x][i][j - i] << std::endl;
				}
			}
			return 0;
		}
		else if (algo == std::string("RintCwithDFTInterval")) {
			typedef WideComplexNumber<IntervalVar> Comp;
			options.allow_fft = false;
			const auto ans1 = ComputeRintW1Dim<Comp>(options);
			const auto ans2 = RintW1DimToBppm(ans1.first, ans1.second);
			for (int x = 0; x < ans2.first.size(); ++x) {
				std::cout << x << " " << ans2.first[x].lower() << " " << ans2.first[x].upper() << std::endl;
				if (ans2.first[x] == 0.0)continue;
				for (int i = 1; i <= n; ++i)for (int j = i + 1; j <= n && j - i <= options.max_span; ++j) {
					std::cout << x << " " << i << " " << j << " " << ans2.second[x][i][j - i].lower() << " " << ans2.second[x][i][j - i].upper() << std::endl;
				}
			}
			return 0;
		}
		else if (algo == std::string("RintCwithFFTInterval")) {
			typedef WideComplexNumber<IntervalVar> Comp;
			options.allow_fft = true;
			const auto ans1 = ComputeRintW1Dim<Comp>(options);
			const auto ans2 = RintW1DimToBppm(ans1.first, ans1.second);
			for (int x = 0; x < ans2.first.size(); ++x) {
				std::cout << x << " " << ans2.first[x].lower() << " " << ans2.first[x].upper() << std::endl;
				if (ans2.first[x] == 0.0)continue;
				for (int i = 1; i <= n; ++i)for (int j = i + 1; j <= n && j - i <= options.max_span; ++j) {
					std::cout << x << " " << i << " " << j << " " << ans2.second[x][i][j - i].lower() << " " << ans2.second[x][i][j - i].upper() << std::endl;
				}
			}
			return 0;
		}
		else if (algo == std::string("RintCwithoutFourierInterval")) {
			typedef WideComplexNumber<IntervalVar> Comp;
			options.allow_fft = false;
			const auto ans1 = ComputeRintW1DimNonFourier<Comp>(options);
			const auto ans2 = RintW1DimToBppm(ans1.first, ans1.second);
			for (int x = 0; x < ans2.first.size(); ++x) {
				std::cout << x << " " << ans2.first[x].lower() << " " << ans2.first[x].upper() << std::endl;
				if (ans2.first[x] == 0.0)continue;
				for (int i = 1; i <= n; ++i)for (int j = i + 1; j <= n && j - i <= options.max_span; ++j) {
					std::cout << x << " " << i << " " << j << " " << ans2.second[x][i][j - i].lower() << " " << ans2.second[x][i][j - i].upper() << std::endl;
				}
			}
			return 0;
		}
	}

	std::cout << "error" << std::endl;
	return 0;

}

}



int main(int argc, char *argv[]) {

	rintdwr::main_(argc, argv);



	return 0;
}
