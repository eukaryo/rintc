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

int verification(const std::string& structure, const std::string& sequence, const int max_span, const int max_loop) {

	for (const char c : sequence) {
		if (!(c == 'A' || c == 'U' || c == 'G' || c == 'C')) {
			std::cerr << "Error: The RNA sequence must be uppercase." << std::endl;
			return 1;
		}
	}
	for (const char c : structure) {
		if (!(c == '(' || c == '.' || c == ')')) {
			std::cerr << "Error: The RNA structure must consist of '(', '.' and ')' only." << std::endl;
			return 1;
		}
	}
	if (!(structure.size() == sequence.size())) {
		std::cerr << "Error: The RNA sequence and structure must be the same length." << std::endl;
		std::cerr << "Your sequence length  = " + std::to_string(sequence.size()) << std::endl;
		std::cerr << "Your structure length = " + std::to_string(structure.size()) << std::endl;
		return 1;
	}
	if (!(1 <= max_span)) {
		std::cerr << "Error: The value of max-span constraint must be a positive integer." << std::endl;
		return 1;
	}
	if (!(max_span <= int(structure.size()))) {
		std::cerr << "Error: The value of max-span constraint must be less than or equal to the length of the RNA sequence." << std::endl;
		std::cerr << "Your sequence length  = " + std::to_string(sequence.size()) << std::endl;
		std::cerr << "Your constraint value = " + std::to_string(max_span) << std::endl;
		return 1;
	}

	{
		const int n = int(sequence.size());
		const std::string bp = "AU UA GC CG GU UG";
		std::string query = "XX";
		std::vector<std::vector<int>>ans(n + 1, std::vector<int>(n + 1, 0));
		std::stack<int> bp_pos;
		for (int i = 1; i <= n; ++i) {
			switch (structure[i - 1]) {
			case '(':
				bp_pos.push(i);
				break;
			case ')':
				if (!(bp_pos.size() >= 1)) {
					std::cerr << "Error: The RNA structure is invalid." << std::endl;
					std::cerr << "')' of position " + std::to_string(i) + " cannot form a base pair." << std::endl;
					return 1;
				}
				if (!(TURN < (i - bp_pos.top()) && (i - bp_pos.top()) <= max_span)) {
					std::cerr << "Error: The RNA structure contains the base pair whose length is longer than the max-span constraint." << std::endl;
					return 1;
				}

				query[0] = sequence[bp_pos.top() - 1];
				query[1] = sequence[i - 1];
				if (!(bp.find(query) != std::string::npos)) {
					std::cerr << "Error: The RNA structure contains an illegal base pair." << std::endl;
					std::cerr << "Position " + std::to_string(bp_pos.top()) + " (the base is '" + query.substr(0, 1) + "') and" << std::endl;
					std::cerr << "position " + std::to_string(i) + " (the base is '" + query.substr(1, 1) + "') " << std::endl;
					return 1;
				}


				ans[bp_pos.top()][i] = 1;
				bp_pos.pop();
				break;
			case '.':
				break;
			default:
				assert(0);
				break;
			}
		}
		if (!(bp_pos.size() == 0)) {
			std::cerr << "Error: The RNA structure is invalid." << std::endl;
			std::cerr << "'(' is more than ')'." << std::endl;
			return 1;
		}
		if (!(ComputeMaxLoop(structure) <= max_loop)) {
			std::cerr << "Error: The RNA structure contains the base pair whose length is longer than the max-loop constraint." << std::endl;
			return 1;
		}
	}
	return 0;
}

int main_(int argc, char *argv[]) {

//#ifdef _WIN64 
//	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
//#endif


	if (argc == 2) {
		if (std::string(argv[1]) == std::string("test")) {
			TestAll();
			return 0;
		}
		if (std::string(argv[1]) == std::string("ComputationTimeExperimentW100")) {
			ComputationTimeExperiment2();
			return 0;
		}
		if (std::string(argv[1]) == std::string("ComputationTimeExperimentWN")) {
			ComputationTimeExperiment3();
			return 0;
		}
		if (std::string(argv[1]) == std::string("HeatResistanceExperiment")) {
			HeatResistanceExperiment("stability_dataset/ecoli_4v9d_ndb_aa.txt", 0.5);
			HeatResistanceExperiment("stability_dataset/thermophilus_4v51_ndb_aa.txt", 0.5);
			HeatResistanceExperimentPK("stability_dataset/ecoli_4v9d_ndb_aa.txt", "stability_dataset/ecoli_4v9d_ndb_aa_structure.txt", 0.5);
			HeatResistanceExperimentPK("stability_dataset/thermophilus_4v51_ndb_aa.txt", "stability_dataset/thermophilus_4v51_ndb_aa.txt", 0.5);
			return 0;
		}
		if (std::string(argv[1]) == std::string("AccuracyExperiment151")) {
			AccuracyExperiment2();
			return 0;
		}
		if (std::string(argv[1]) == std::string("RrnaAccuracyExperiment")) {
			RrnaAccuracyExperiment("stability_dataset/ecoli_4v9d_ndb_aa.txt");
			RrnaAccuracyExperiment("stability_dataset/thermophilus_4v51_ndb_aa.txt");
			return 0;
		}
	}
	if (argc == 4) {
		if (std::string(argv[1]) == std::string("ComputationTimeExperimentWvar")) {
			const int i = std::stoi(std::string(argv[2]));
			const int j = std::stoi(std::string(argv[3]));
			omp_set_num_threads(std::min<int>(j, omp_get_max_threads()));
			switch (i) {
			case 6: ComputationTimeExperiment5(50); return 0;
			case 5: ComputationTimeExperiment5(100); return 0;
			case 4: ComputationTimeExperiment5(200); return 0;
			case 3: ComputationTimeExperiment5(300); return 0;
			case 2: ComputationTimeExperiment5(400); return 0;
			case 1: ComputationTimeExperiment5(10000); return 0;
			}
			return 0;
		}
	}
	if (argc == 5) {

		const std::string sequence = std::string(argv[1]);
		const std::string structure = std::string(argv[2]);
		const int W = std::stoi(std::string(argv[3]));
		const std::string algo = std::string(argv[4]);
		const int n = sequence.length();
		const int max_loop = n < 30 ? n : 30;

		if (verification(structure, sequence, W, max_loop)) {
			return 1;
		}

		RintX1DOptions options;
		options.temperature = 37.0;
		options.param_file_name = std::string("Turner2004");
		options.max_span = W;
		options.max_loop = max_loop;
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
		std::cout << "Error: invalid algo." << std::endl;
		return 1;
	}

	std::cout << "Error" << std::endl;
	return 1;


}

}



int main(int argc, char *argv[]) {

	return rintdwr::main_(argc, argv);

}
