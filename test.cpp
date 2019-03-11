/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#include"test.h"

#include"complex_number.h"
#include"interval_type.h"
#include"parameter.h"
#include"misc.h"
#include"rintx.h"
#include"mccaskill_1990.h"
#include"sample_mccaskill.h"
#include"centroid_fold.h"

namespace rintdwr {

void TestMaxHamming(const int num) {

	std::mt19937_64 rnd(123456);
	std::uniform_int_distribution<int> length(15, 22);
	std::uniform_int_distribution<int> span(0, 13);

	for (int iteration = 0; iteration < num; ++iteration) {
		std::cout << iteration << std::endl;

		//ComputeMaxHammingDistanceが正しい解を返すことを確認する。
		const int n = length(rnd);
		const int max_span = (iteration % 2) ? n - span(rnd) : n;
		const int max_loop = ((iteration / 2) % 2) ? 13 - span(rnd) : n;
		const std::string sequence = MakeRandomRnaSeq(n, rnd());
		const std::vector<std::string>structures = EnumerateStructures(sequence, max_span, max_loop);
		std::uniform_int_distribution<int> length(0, int(structures.size()) - 1);
		const std::string reference_structure = structures[length(rnd)];
		const auto base_pairing_matrix = VerifyAndParseStructure(reference_structure, sequence, max_span, max_loop);
		const int max_dim = ComputeMaxHammingDistance(sequence, base_pairing_matrix, max_span, max_loop);

		int maxd = 0;
		for (const auto s : structures)maxd = std::max<int>(maxd, ComputeHammingDistance(reference_structure, s));
		assert(maxd == max_dim);
	}
}

void TestSampling(const int num) {

	std::mt19937_64 rnd(123456);
	std::uniform_int_distribution<int> length(15, 22);
	std::uniform_int_distribution<int> span(0, 13);

	std::uniform_int_distribution<int> longlength(30, 50);
	std::uniform_int_distribution<int> longspan(1, 40);

	for (int iteration = 0; iteration < num; ++iteration) {
		std::cout << iteration << std::endl;

		{
			//すごく短いRNAの配列に対して、
			//Brute Forceで全列挙した構造の数と内側分配関数の値が等しいことを確かめる。
			const int n = length(rnd);
			const int max_span = (iteration % 2) ? n - span(rnd) : n;
			const int max_loop = ((iteration / 2) % 2) ? 13 - span(rnd) : n;
			const std::string sequence = MakeRandomRnaSeq(n, rnd());
			const std::vector<std::string>structures = EnumerateStructures(sequence, max_span, max_loop);

			const auto ans1 = SampleMcCaskillUniform(sequence, 1, max_span, max_loop, true, int(rnd()));
			assert(ans1.second == Bigint(structures.size()));

			const auto ans2 = SampleMcCaskillEnergyAware(sequence, "Turner2004", 37.0, 1, max_span, max_loop, true, int(rnd()));
			assert(abs(ans2.second.log_scale - log(SimpleMcCaskill<double>(sequence, "Turner2004", 37.0, max_span, max_loop).second)) < 0.00001);

		}

		{
			//長めのRNA配列に対して、単なるMcCaskill型DPとSamplingで塩基対確率行列がほぼ一致することを確かめる。
			const int n = longlength(rnd);
			const int max_span = (iteration % 2) ? std::min(n, longspan(rnd)) : n;
			const int max_loop = ((iteration / 2) % 2) ? (longspan(rnd) % 30) : 30;
			const std::string sequence = MakeRandomRnaSeq(n, rnd());
			const auto ans1 = SampleMcCaskillUniform(sequence, 1, max_span, max_loop, true, int(rnd()));
			const auto ans2 = SampleMcCaskillEnergyAware(sequence, "Turner2004", 37.0, 1, max_span, max_loop, true, int(rnd()));
		}
	}
}

void TestSimpleMcCaskillWide(const int num) {

	std::mt19937_64 rnd(123456);

	std::uniform_int_distribution<int> longlength(40, 70);
	std::uniform_int_distribution<int> longspan(1, 60);

	RintX1DOptions options;
	options.temperature = 37.0;
	options.param_file_name = std::string("Turner2004");

	for (int iteration = 0; iteration < num; ++iteration) {
		std::cout << iteration << std::endl;

		{
			//長めのRNA配列に対して、WideとUsualで内側分配関数が一致することを確かめる。
			const int n = longlength(rnd);
			options.max_span = (iteration % 2) ? std::min(n, longspan(rnd)) : n;
			options.max_loop = ((iteration / 2) % 2) ? (longspan(rnd) % 30) : 30;
			options.allow_fft = ((iteration / 4) % 2) == 0;
			options.sequence = MakeRandomRnaSeq(n, rnd());
			options.reference_structure1 = SampleMcCaskillUniform(options.sequence, 1, options.max_span, options.max_loop, false, int(rnd())).first[0];
			options.S1 = VerifyAndParseStructure(options.reference_structure1, options.sequence, options.max_span, options.max_loop);
			options.max_dim1 = ComputeMaxHammingDistance(options.sequence, options.S1, options.max_span, options.max_loop);
			const auto ans1 = SimpleMcCaskill<Floating>(options.sequence, options.param_file_name, options.temperature, options.max_span, options.max_loop);
			const auto ans2 = SimpleMcCaskillWide(options.sequence, options.param_file_name, options.temperature, options.max_span, options.max_loop);

			for (int i = 1; i <= n; ++i)for (int j = i + 1; j <= n && j - i <= options.max_span; ++j) {
				const Floating diff = abs(ans1.first[i][j - i] - ans2.first[i][j - i]);
				assert(diff < Floating(0.0000000001));
			}
			assert(abs(log(ans1.second) - ans2.second) < Floating(0.0000000001));
		}
	}
}

void TestRintD1Dim(const int num) {

	typedef WideComplexNumber<Floating> Comp;

	std::mt19937_64 rnd(123456);
	std::uniform_int_distribution<int> length(15, 22);
	std::uniform_int_distribution<int> span(0, 13);

	std::uniform_int_distribution<int> longlength(40, 70);
	std::uniform_int_distribution<int> longspan(1, 60);

	RintX1DOptions options;
	options.temperature = 37.0;
	options.param_file_name = std::string("Turner2004");

	for (int iteration = 0; iteration < num; ++iteration) {
		std::cout << iteration << std::endl;

		{
			//すごく短いRNA配列に対して、Brute ForceとDP+DFTで解が一致することを確かめる。
			const int n = length(rnd);
			options.max_span = (iteration % 2) ? n - span(rnd) : n;
			options.max_loop = ((iteration / 2) % 2) ? 13 - span(rnd) : n;
			options.allow_fft = ((iteration / 4) % 2) == 0;
			options.sequence = MakeRandomRnaSeq(n, rnd());
			const std::vector<std::string>structures = EnumerateStructures(options.sequence, options.max_span, options.max_loop);
			std::uniform_int_distribution<int> length(0, int(structures.size()) - 1);
			options.reference_structure1 = structures[length(rnd)];
			options.S1 = VerifyAndParseStructure(options.reference_structure1, options.sequence, options.max_span, options.max_loop);
			options.max_dim1 = ComputeMaxHammingDistance(options.sequence, options.S1, options.max_span, options.max_loop);
			const auto ans1 = BruteForceRintD1Dim(options);
			const auto ans2 = ComputeRintD1Dim<Comp>(options);

			double maxnum = 0.0;
			double maxdiff = 0.0;
			for (int d = 0; d <= options.max_dim1; ++d) {
				maxnum = std::max<double>(maxnum, abs(ans1[d]));
				maxdiff = std::max<double>(maxdiff, abs(ans1[d] - double(ans2[d].ToUsualComp().real)));
			}
			if (maxdiff >= maxnum * 0.0000000001) {
				options.Dump();
				assert(0);
			}
		}

		{
			//長めのRNA配列に対して、単なるMcCaskill型DPとMoriDP+DFTで内側分配関数が一致することを確かめる。
			const int n = longlength(rnd);
			options.max_span = (iteration % 2) ? std::min(n, longspan(rnd)) : n;
			options.max_loop = ((iteration / 2) % 2) ? (longspan(rnd) % 30) : 30;
			options.allow_fft = ((iteration / 4) % 2) == 0;
			options.sequence = MakeRandomRnaSeq(n, rnd());
			options.reference_structure1 = SampleMcCaskillUniform(options.sequence, 1, options.max_span, options.max_loop, false, int(rnd())).first[0];
			options.S1 = VerifyAndParseStructure(options.reference_structure1, options.sequence, options.max_span, options.max_loop);
			options.max_dim1 = ComputeMaxHammingDistance(options.sequence, options.S1, options.max_span, options.max_loop);
			const auto ans1 = SimpleMcCaskill<Floating>(options.sequence, options.param_file_name, options.temperature, options.max_span, options.max_loop);
			const auto ans2 = ComputeRintD1Dim<Comp>(options);

			double max_scale = 0.0;
			for (const auto c : ans2) max_scale = std::max<double>(max_scale, double(c.log_scale));

			double mori_sum = 0.0;
			for (const auto c : ans2) {
				if (c.IsZero())continue;
				assert(abs(c.real - 1.0) * exp(c.log_scale - max_scale) < 0.0000000001);
				assert(abs(c.imag) * exp(c.log_scale - max_scale) < 0.0000000001);
				mori_sum += exp(double(c.log_scale));
			}
			assert(abs(mori_sum - double(ans1.second)) < std::max<double>(mori_sum, double(ans1.second)) * 0.0000000001);
		}
	}
}

void TestRintD2Dim(const int num) {

	typedef WideComplexNumber<Floating> Comp;

	std::mt19937_64 rnd(123456);
	std::uniform_int_distribution<int> length(15, 22);
	std::uniform_int_distribution<int> span(0, 13);

	std::uniform_int_distribution<int> longlength(30, 60);
	std::uniform_int_distribution<int> longspan(1, 40);

	RintX2DOptions options;
	options.temperature = 37.0;
	options.param_file_name = std::string("Turner2004");

	for (int iteration = 0; iteration < num; ++iteration) {
		std::cout << iteration << std::endl;

		{
			//すごく短いRNA配列に対して、Brute ForceとDP+DFTで解が一致することを確かめる。
			const int n = length(rnd);
			options.max_span = (iteration % 2) ? n - span(rnd) : n;
			options.max_loop = ((iteration / 2) % 2) ? 13 - span(rnd) : n;
			options.allow_fft = ((iteration / 4) % 2) == 0;
			options.sequence = MakeRandomRnaSeq(n, rnd());
			const std::vector<std::string>structures = EnumerateStructures(options.sequence, options.max_span,options.max_loop);
			std::uniform_int_distribution<int> length(0, int(structures.size()) - 1);
			options.reference_structure1 = structures[length(rnd)];
			options.reference_structure2 = structures[length(rnd)];
			options.S1 = VerifyAndParseStructure(options.reference_structure1, options.sequence, options.max_span, options.max_loop);
			options.S2 = VerifyAndParseStructure(options.reference_structure2, options.sequence, options.max_span, options.max_loop);
			options.max_dim1 = ComputeMaxHammingDistance(options.sequence, options.S1, options.max_span, options.max_loop);
			options.max_dim2 = ComputeMaxHammingDistance(options.sequence, options.S2, options.max_span, options.max_loop);
			const auto ans1 = BruteForceRintD2Dim(options);
			const auto ans2 = ComputeRintD2Dim<Comp>(options);

			double maxnum = 0.0;
			double maxdiff = 0.0;
			for (int d1 = 0; d1 <= options.max_dim1; ++d1)for (int d2 = 0; d2 <= options.max_dim2; ++d2) {
				maxnum = std::max<double>(maxnum, abs(ans1[d1][d2]));
				maxdiff = std::max<double>(maxdiff, double(abs(ans1[d1][d2] - ans2[d1][d2].ToUsualComp().real)));
			}
			if (maxdiff >= maxnum * 0.0000000001) {
				options.Dump();
				assert(0);
			}
		}

		{
			//長めのRNA配列に対して、単なるMcCaskill型DPとMoriDP+DFTで内側分配関数が一致することを確かめる。
			const int n = longlength(rnd);
			options.max_span = (iteration % 2) ? std::min(n, longspan(rnd)) : n;
			options.max_loop = ((iteration / 2) % 2) ? (longspan(rnd) % 30) : 30;
			options.allow_fft = ((iteration / 4) % 2) == 0;
			options.sequence = MakeRandomRnaSeq(n, rnd());
			options.reference_structure1 = SampleMcCaskillUniform(options.sequence, 1, options.max_span, options.max_loop, false, int(rnd())).first[0];
			options.reference_structure2 = SampleMcCaskillUniform(options.sequence, 1, options.max_span, options.max_loop, false, int(rnd())).first[0];
			options.S1 = VerifyAndParseStructure(options.reference_structure1, options.sequence, options.max_span, options.max_loop);
			options.S2 = VerifyAndParseStructure(options.reference_structure2, options.sequence, options.max_span, options.max_loop);
			options.max_dim1 = ComputeMaxHammingDistance(options.sequence, options.S1, options.max_span, options.max_loop);
			options.max_dim2 = ComputeMaxHammingDistance(options.sequence, options.S2, options.max_span, options.max_loop);
			const auto ans1 = SimpleMcCaskill<Floating>(options.sequence, options.param_file_name, options.temperature, options.max_span, options.max_loop);
			const auto ans2 = ComputeRintD2Dim<Comp>(options);

			double max_scale = 0.0;
			for (int i = 0; i <= options.max_dim1; ++i)for (int j = 0; j <= options.max_dim2; ++j) {
				max_scale = std::max<double>(max_scale, double(ans2[i][j].log_scale));
			}

			double mori_sum = 0.0;
			for (int i = 0; i <= options.max_dim1; ++i)for (int j = 0; j <= options.max_dim2; ++j) {
				if (ans2[i][j].IsZero())continue;
				assert(abs(ans2[i][j].real - 1.0) * exp(ans2[i][j].log_scale - max_scale) < 0.0000000001);
				assert(abs(ans2[i][j].imag) * exp(ans2[i][j].log_scale - max_scale) < 0.0000000001);
				mori_sum += exp(double(ans2[i][j].log_scale));
			}
			assert(abs(mori_sum - double(ans1.second)) < std::max<double>(mori_sum, double(ans1.second)) * 0.0000000001);
		}
	}

}

void TestRintW1Dim(const int num) {

	typedef WideComplexNumber<Floating> Comp;

	std::mt19937_64 rnd(123456);
	std::uniform_int_distribution<int> length(15, 22);
	std::uniform_int_distribution<int> span(0, 13);

	std::uniform_int_distribution<int> longlength(40, 70);
	std::uniform_int_distribution<int> longspan(1, 60);

	RintX1DOptions options;
	options.temperature = 37.0;
	options.param_file_name = std::string("Turner2004");

	for (int iteration = 0; iteration < num; ++iteration) {
		std::cout << iteration << std::endl;

		{
			//すごく短いRNA配列に対して、Brute ForceとDP+DFTで解が一致することを確かめる。

			const int n = length(rnd);
			options.max_span = (iteration % 2) ? n - span(rnd) : n;
			options.max_loop = ((iteration / 2) % 2) ? 13 - span(rnd) : n;
			options.allow_fft = ((iteration / 4) % 2) == 0;
			options.sequence = MakeRandomRnaSeq(n, rnd());
			const std::vector<std::string>structures = EnumerateStructures(options.sequence, options.max_span, options.max_loop);
			std::uniform_int_distribution<int> length(0, int(structures.size()) - 1);
			options.reference_structure1 = structures[length(rnd)];
			options.S1 = VerifyAndParseStructure(options.reference_structure1, options.sequence, options.max_span, options.max_loop);
			options.max_dim1 = ComputeMaxHammingDistance(options.sequence, options.S1, options.max_span, options.max_loop);
			const auto ans1 = BruteForceRintW1Dim(options);
			const auto ans2 = ComputeRintW1Dim<Comp>(options);

			double maxnum = 0.0;
			double maxdiff = 0.0;
			for (int d = 0; d <= options.max_dim1; ++d) {
				maxnum = std::max<double>(maxnum, abs(ans1.first[d]));
				maxdiff = std::max<double>(maxdiff, double(abs(ans1.first[d] - ans2.first[d].ToUsualComp().real)));
			}
			for (int d = 0; d <= options.max_dim1; ++d) {
				for (int i = 1; i <= n; i++)for (int j = i; j <= n && j - i <= options.max_span; j++) {
					maxnum = std::max<double>(maxnum, abs(ans1.second[d][i][j - i]));
					maxdiff = std::max<double>(maxdiff, double(abs(ans1.second[d][i][j - i] - ans2.second[d][i][j - i].ToUsualComp().real)));
				}
			}
			if (maxdiff >= maxnum * 0.0000000001) {
				options.Dump();
				assert(0);
			}
		}

		{
			//長めのRNA配列に対して、単なるMcCaskill型DPとHagioDP+DFTで内側分配関数と塩基対確率行列が一致することを確かめる。
			const int n = longlength(rnd);
			options.max_span = (iteration % 2) ? std::min(n, longspan(rnd)) : n;
			options.max_loop = ((iteration / 2) % 2) ? (longspan(rnd) % 30) : 30;
			options.allow_fft = ((iteration / 4) % 2) == 0;
			options.sequence = MakeRandomRnaSeq(n, rnd());
			options.reference_structure1 = SampleMcCaskillUniform(options.sequence, 1, options.max_span, options.max_loop, false, int(rnd())).first[0];
			options.S1 = VerifyAndParseStructure(options.reference_structure1, options.sequence, options.max_span, options.max_loop);
			options.max_dim1 = ComputeMaxHammingDistance(options.sequence, options.S1, options.max_span, options.max_loop);
			const auto ans1 = SimpleMcCaskill<Floating>(options.sequence, options.param_file_name, options.temperature, options.max_span, options.max_loop);
			const auto ans2 = ComputeRintW1Dim<Comp>(options);
			const auto ans3 = RintW1DimToBppm(ans2.first, ans2.second);

			double max_scale = 0.0;
			for (const auto c : ans2.first) max_scale = std::max<double>(max_scale, double(c.log_scale));

			double rintd_sum = 0.0;
			for (const auto c : ans2.first) {
				if (c.IsZero())continue;
				assert(abs(c.real - 1.0) * exp(c.log_scale - max_scale) < 0.0000000001);
				assert(abs(c.imag) * exp(c.log_scale - max_scale) < 0.0000000001);
				rintd_sum += exp(double(c.log_scale));
			}
			assert(abs(rintd_sum - double(ans1.second)) < std::max<double>(rintd_sum, double(ans1.second)) * 0.0000000001);

			for (int i = 1; i <= n; ++i)for (int j = i + 1; j <= n && j - i <= options.max_span; ++j) {
				double rintw_prob = 0.0;
				for (int x = 0; x <= options.max_dim1; ++x)rintw_prob += double(ans3.first[x] * ans3.second[x][i][j - i]);
				assert(abs(ans1.first[i][j - i] - rintw_prob) < 0.01);
			}
		}
	}
}

void TestRintW2Dim(const int num) {

	typedef WideComplexNumber<Floating> Comp;

	std::mt19937_64 rnd(123456);
	std::uniform_int_distribution<int> length(15, 22);
	std::uniform_int_distribution<int> span(0, 13);

	std::uniform_int_distribution<int> longlength(30, 60);
	std::uniform_int_distribution<int> longspan(1, 40);

	RintX2DOptions options;
	options.temperature = 37.0;
	options.param_file_name = std::string("Turner2004");

	for (int iteration = 0; iteration < num; ++iteration) {
		std::cout << iteration << std::endl;

		{
			//すごく短いRNA配列に対して、Brute ForceとDP+DFTで解が一致することを確かめる。

			const int n = length(rnd);
			options.max_span = (iteration % 2) ? n - span(rnd) : n;
			options.max_loop = ((iteration / 2) % 2) ? 13 - span(rnd) : n;
			options.allow_fft = ((iteration / 4) % 2) == 0;
			options.sequence = MakeRandomRnaSeq(n, rnd());
			const std::vector<std::string>structures = EnumerateStructures(options.sequence, options.max_span, options.max_loop);
			std::uniform_int_distribution<int> length(0, int(structures.size()) - 1);
			options.reference_structure1 = structures[length(rnd)];
			options.reference_structure2 = structures[length(rnd)];
			options.S1 = VerifyAndParseStructure(options.reference_structure1, options.sequence, options.max_span, options.max_loop);
			options.S2 = VerifyAndParseStructure(options.reference_structure2, options.sequence, options.max_span, options.max_loop);
			options.max_dim1 = ComputeMaxHammingDistance(options.sequence, options.S1, options.max_span, options.max_loop);
			options.max_dim2 = ComputeMaxHammingDistance(options.sequence, options.S2, options.max_span, options.max_loop);
			const auto ans1 = BruteForceRintW2Dim(options);
			const auto ans2 = ComputeRintW2Dim<Comp>(options);

			double maxnum = 0.0;
			double maxdiff = 0.0;
			for (int d1 = 0; d1 <= options.max_dim1; ++d1)for (int d2 = 0; d2 <= options.max_dim2; ++d2) {
				maxnum = std::max<double>(maxnum, abs(ans1.first[d1][d2]));
				maxdiff = std::max<double>(maxdiff, double(abs(ans1.first[d1][d2] - ans2.first[d1][d2].ToUsualComp().real)));
			}
			for (int d1 = 0; d1 <= options.max_dim1; ++d1)for (int d2 = 0; d2 <= options.max_dim2; ++d2) {
				for (int i = 1; i <= n; i++)for (int j = i; j <= n && j - i <= options.max_span; j++) {
					maxnum = std::max<double>(maxnum, abs(ans1.second[d1][d2][i][j - i]));
					maxdiff = std::max<double>(maxdiff, double(abs(ans1.second[d1][d2][i][j - i] - ans2.second[d1][d2][i][j - i].ToUsualComp().real)));
				}
			}
			if (maxdiff >= maxnum * 0.0000000001) {
				options.Dump();
				assert(0);
			}
		}

		{
			//長めのRNA配列に対して、単なるMcCaskill型DPとHagioDP+DFTで内側分配関数と塩基対確率行列が一致することを確かめる。
			const int n = longlength(rnd);
			options.max_span = (iteration % 2) ? std::min(n, longspan(rnd)) : n;
			options.max_loop = ((iteration / 2) % 2) ? (longspan(rnd) % 30) : 30;
			options.allow_fft = ((iteration / 4) % 2) == 0;
			options.sequence = MakeRandomRnaSeq(n, rnd());
			options.reference_structure1 = SampleMcCaskillUniform(options.sequence, 1, options.max_span, options.max_loop, false, int(rnd())).first[0];
			options.reference_structure2 = SampleMcCaskillUniform(options.sequence, 1, options.max_span, options.max_loop, false, int(rnd())).first[0];
			options.S1 = VerifyAndParseStructure(options.reference_structure1, options.sequence, options.max_span, options.max_loop);
			options.S2 = VerifyAndParseStructure(options.reference_structure2, options.sequence, options.max_span, options.max_loop);
			options.max_dim1 = ComputeMaxHammingDistance(options.sequence, options.S1, options.max_span, options.max_loop);
			options.max_dim2 = ComputeMaxHammingDistance(options.sequence, options.S2, options.max_span, options.max_loop);
			const auto ans1 = SimpleMcCaskill<Floating>(options.sequence, options.param_file_name, options.temperature, options.max_span, options.max_loop);
			const auto ans2 = ComputeRintW2Dim<Comp>(options);
			const auto ans3 = RintW2DimToBppm(ans2.first, ans2.second);

			double max_scale = 0.0;
			for (int d1 = 0; d1 <= options.max_dim1; ++d1)for (int d2 = 0; d2 <= options.max_dim2; ++d2) {
				max_scale = std::max<double>(max_scale, double(ans2.first[d1][d2].log_scale));
			}

			double rintd_sum = 0.0;
			for (int d1 = 0; d1 <= options.max_dim1; ++d1)for (int d2 = 0; d2 <= options.max_dim2; ++d2) {
				if (ans2.first[d1][d2].IsZero())continue;
				assert(abs(ans2.first[d1][d2].real - 1.0) * exp(ans2.first[d1][d2].log_scale - max_scale) < 0.0000000001);
				assert(abs(ans2.first[d1][d2].imag) * exp(ans2.first[d1][d2].log_scale - max_scale) < 0.0000000001);
				rintd_sum += exp(double(ans2.first[d1][d2].log_scale));
			}
			assert(abs(rintd_sum - double(ans1.second)) < std::max<double>(rintd_sum, double(ans1.second)) * 0.0000000001);

			for (int i = 1; i <= n; ++i)for (int j = i + 1; j <= n && j - i <= options.max_span; ++j) {
				double rintw_prob = 0.0;
				for (int d1 = 0; d1 <= options.max_dim1; ++d1)for (int d2 = 0; d2 <= options.max_dim2 && j - i <= options.max_span; ++d2) {
					rintw_prob += double(ans3.first[d1][d2] * ans3.second[d1][d2][i][j - i]);
				}
				assert(abs(ans1.first[i][j - i] - rintw_prob) < 0.01);
			}
		}
	}
}

void TestCentroidFold(const int num) {

	typedef WideComplexNumber<Floating> Comp;

	std::mt19937_64 rnd(123457);
	std::uniform_int_distribution<int> length(15, 20);
	std::uniform_int_distribution<int> span(0, 13);

	std::uniform_int_distribution<int> longlength(30, 65);
	std::uniform_int_distribution<int> longspan(1, 50);

	RintX1DOptions options;
	options.temperature = 37.0;
	options.param_file_name = std::string("Turner2004");

	std::uniform_real_distribution<double>gamma_rnd(0.1, 3.0);

	for (int iteration = 0; iteration < num; ++iteration) {
		std::cout << iteration << std::endl;

		const Floating gamma = Floating(gamma_rnd(rnd));

		{
			//短いRNA配列に対して、Brute Forceで求めたMEG構造と合っていることを確認する。

			const int n = length(rnd);
			options.max_span = (iteration % 2) ? n - span(rnd) : n;
			options.max_loop = ((iteration / 2) % 2) ? 13 - span(rnd) : n;
			options.allow_fft = ((iteration / 4) % 2) == 0;
			options.sequence = MakeRandomRnaSeq(n, rnd());
			const std::vector<std::string>structures = EnumerateStructures(options.sequence, options.max_span, options.max_loop);
			std::uniform_int_distribution<int> length(0, int(structures.size()) - 1);
			options.reference_structure1 = structures[length(rnd)];
			options.S1 = VerifyAndParseStructure(options.reference_structure1, options.sequence, options.max_span, options.max_loop);
			options.max_dim1 = ComputeMaxHammingDistance(options.sequence, options.S1, options.max_span, options.max_loop);
			const auto ans1 = SimpleMcCaskill<Floating>(options.sequence, options.param_file_name, options.temperature, options.max_span, options.max_loop);
			const auto ans2 = ComputeRintW1Dim<Comp>(options);

			auto ans = RintW1DimToBppm(ans2.first, ans2.second).second;
			ans.push_back(ans1.first);

			for (const auto c : ans) {

				const auto fold1 = GetCentroidFoldMcCaskill(options.sequence, options.max_span, c, gamma, options.max_loop);
				const auto fold2 = BruteForceGetCentroidFoldMcCaskill(options.sequence, options.max_span, c, gamma, options.max_loop);

				const double diff = ComputeGain(fold1, c, gamma, options.max_loop) - ComputeGain(fold2.first, c, gamma, options.max_loop);
				const double maxnum = std::max<double>(abs(ComputeGain(fold1, c, gamma, options.max_loop)), abs(ComputeGain(fold2.first, c, gamma, options.max_loop)));
				assert(abs(diff) <= maxnum * 0.00000000001);

				const auto folds1 = GetCentroidFoldForEachHammingDistance(options.sequence, options.S1, options.max_dim1, options.max_span, c, gamma, options.max_loop);
				const auto folds2 = BruteForceGetCentroidFoldForEachHammingDistance(options.sequence, options.reference_structure1, options.max_dim1, options.max_span, c, gamma, options.max_loop);

				for (int i = 0; i <= options.max_dim1; ++i) {
					const double diff = ComputeGain(folds1[i], c, gamma, options.max_loop) - ComputeGain(folds2[i].first, c, gamma, options.max_loop);
					const double maxnum = std::max<double>(abs(ComputeGain(folds1[i], c, gamma, options.max_loop)), abs(ComputeGain(folds2[i].first, c, gamma, options.max_loop)));
					assert(abs(diff) <= maxnum * 0.00000000001);
				}
			}
		}

		{
			//長めのRNA配列に対して、答えの形式が正しいことと、ハミング距離ごとの解に真のMEG構造が含まれることを確認する。

			const int n = longlength(rnd);
			options.max_span = (iteration % 2) ? std::min(n, longspan(rnd)) : n;
			options.max_loop = ((iteration / 2) % 2) ? (longspan(rnd) % 30) : 30;
			options.allow_fft = ((iteration / 4) % 2) == 0;
			options.sequence = MakeRandomRnaSeq(n, rnd());
			options.reference_structure1 = SampleMcCaskillUniform(options.sequence, 1, options.max_span, options.max_loop, false, int(rnd())).first[0];
			options.S1 = VerifyAndParseStructure(options.reference_structure1, options.sequence, options.max_span, options.max_loop);
			options.max_dim1 = ComputeMaxHammingDistance(options.sequence, options.S1, options.max_span, options.max_loop);
			const auto ans1 = SimpleMcCaskill<Floating>(options.sequence, options.param_file_name, options.temperature, options.max_span, options.max_loop);
			const auto ans2 = ComputeRintW1Dim<Comp>(options);

			auto ans = RintW1DimToBppm(ans2.first, ans2.second).second;
			ans.push_back(ans1.first);

			for (const auto c : ans) {

				const auto meg = GetCentroidFoldMcCaskill(options.sequence, options.max_span, c, gamma, options.max_loop);
				const int meg_dist = ComputeHammingDistance(options.reference_structure1, meg);

				const auto folds = GetCentroidFoldForEachHammingDistance(options.sequence, options.S1, options.max_dim1, options.max_span, c, gamma, options.max_loop);
				assert(int(folds.size()) == options.max_dim1 + 1);
				for (int i = 0; i <= options.max_dim1; i++) {
					VerifyAndParseStructure(folds[i], options.sequence, options.max_span, options.max_loop);
					assert(ComputeHammingDistance(options.reference_structure1, folds[i]) == i);
					if (meg_dist == i) {
						const double diff = ComputeGain(meg, c, gamma, options.max_loop) - ComputeGain(folds[i], c, gamma, options.max_loop);
						const double maxnum = std::max<double>(abs(ComputeGain(meg, c, gamma, options.max_loop)), abs(ComputeGain(folds[i], c, gamma, options.max_loop)));
						assert(abs(diff) <= maxnum * 0.00000000001);
					}
				}
			}
		}
	}
}

void TestHagioNonFourier(const int num) {

	typedef WideComplexNumber<IntervalVar> Comp;

	std::mt19937_64 rnd(123456);

	std::uniform_int_distribution<int> longlength(20, 40);
	std::uniform_int_distribution<int> longspan(1, 40);

	RintX1DOptions options;
	options.temperature = 37.0;
	options.param_file_name = std::string("Turner2004");

	for (int iteration = 0; iteration < num; ++iteration) {
		std::cout << iteration << std::endl;

		{
			//長めのRNA配列に対して、DFT版とNonFourier版で精度保証区間がoverlapすることを確かめる。
			const int n = longlength(rnd);
			options.max_span = (iteration % 2) ? std::min(n, longspan(rnd)) : n;
			options.max_loop = ((iteration / 2) % 2) ? (longspan(rnd) % 30) : 30;
			options.allow_fft = ((iteration / 4) % 2) == 0;
			options.sequence = MakeRandomRnaSeq(n, rnd());
			options.reference_structure1 = SampleMcCaskillUniform(options.sequence, 1, options.max_span, options.max_loop, false, int(rnd())).first[0];
			options.S1 = VerifyAndParseStructure(options.reference_structure1, options.sequence, options.max_span, options.max_loop);
			options.max_dim1 = ComputeMaxHammingDistance(options.sequence, options.S1, options.max_span, options.max_loop);


			const auto ans3 = ComputeRintW1DimNonFourier<Comp>(options);
			const auto ans4 = ComputeRintW1Dim<Comp>(options);
			for (int x = 0; x <= options.max_dim1; x++) {
				const IntervalVar a1 = ans3.first[x].ToUsualComp().real;
				const IntervalVar a2 = ans4.first[x].ToUsualComp().real;
				if (!overlap(a1, a2)) {
					std::cout << std::setprecision(20) << a1 << std::endl << a2 << std::endl;
				}
				assert(overlap(a1, a2));
			}
			for (int x = 0; x <= options.max_dim1; x++)for (int i = 1; i <= n; ++i)for (int j = i + 1; j <= n && j - i <= options.max_span; ++j) {
				const IntervalVar a1 = ans3.second[x][i][j - i].ToUsualComp().real;
				const IntervalVar a2 = ans4.second[x][i][j - i].ToUsualComp().real;
				if (!overlap(a1, a2)) {
					std::cout << std::setprecision(20) << a1 << std::endl << a2 << std::endl;
				}
				assert(overlap(a1, a2));
			}
		}
	}
}

void TestMaxHammingPK(const int num) {

	std::mt19937_64 rnd(123456);
	std::uniform_int_distribution<int> length(15, 22);
	std::uniform_int_distribution<int> span(0, 13);

	for (int iteration = 0; iteration < num; ++iteration) {
		std::cout << iteration << std::endl;


		//ComputeMaxHammingDistanceが正しい解を返すことを確認する。



		const int n = length(rnd);
		const int max_span = (iteration % 2) ? n - span(rnd) : n;
		const int max_loop = ((iteration / 2) % 2) ? 13 - span(rnd) : n;
		const std::string sequence = MakeRandomRnaSeq(n, rnd());
		const std::vector<std::string>structures = EnumerateStructures(sequence, max_span, max_loop);
		std::uniform_int_distribution<int> length(0, int(structures.size()) - 1);
		const std::vector<std::vector<int>> reference_structure = RandomStructurePK(sequence, rnd());
		const int max_dim = ComputeMaxHammingDistance(sequence, reference_structure, max_span, max_loop);

		int maxd = 0;
		for (const auto s : structures)maxd = std::max<int>(maxd, ComputeHammingDistance(reference_structure, VerifyAndParseStructure(s,sequence,max_span,max_loop)));
		assert(maxd == max_dim);
	}
}

void TestRintD1DimPK(const int num) {

	typedef WideComplexNumber<Floating> Comp;

	std::mt19937_64 rnd(123456);
	std::uniform_int_distribution<int> length(15, 22);
	std::uniform_int_distribution<int> span(0, 13);

	std::uniform_int_distribution<int> longlength(40, 70);
	std::uniform_int_distribution<int> longspan(1, 60);

	RintX1DOptions options;
	options.temperature = 37.0;
	options.param_file_name = std::string("Turner2004");

	for (int iteration = 0; iteration < num; ++iteration) {
		std::cout << iteration << std::endl;

		{
			//すごく短いRNA配列に対して、Brute ForceとDP+DFTで解が一致することを確かめる。
			const int n = length(rnd);
			options.max_span = (iteration % 2) ? n - span(rnd) : n;
			options.max_loop = ((iteration / 2) % 2) ? 13 - span(rnd) : n;
			options.allow_fft = ((iteration / 4) % 2) == 0;
			options.sequence = MakeRandomRnaSeq(n, rnd());
			const std::vector<std::string>structures = EnumerateStructures(options.sequence, options.max_span, options.max_loop);
			std::uniform_int_distribution<int> length(0, int(structures.size()) - 1);
			options.reference_structure1 = std::string("");
			options.S1 = RandomStructurePK(options.sequence, rnd());
			options.max_dim1 = ComputeMaxHammingDistance(options.sequence, options.S1, options.max_span, options.max_loop);
			const auto ans1 = BruteForceRintD1DimPK(options);
			const auto ans2 = ComputeRintD1Dim<Comp>(options);

			double maxnum = 0.0;
			double maxdiff = 0.0;
			for (int d = 0; d <= options.max_dim1; ++d) {
				maxnum = std::max<double>(maxnum, abs(ans1[d]));
				maxdiff = std::max<double>(maxdiff, abs(ans1[d] - double(ans2[d].ToUsualComp().real)));
			}
			if (maxdiff >= maxnum * 0.0000000001) {
				options.Dump();
				assert(0);
			}
		}

		{
			//長めのRNA配列に対して、単なるMcCaskill型DPとMoriDP+DFTで内側分配関数が一致することを確かめる。
			const int n = longlength(rnd);
			options.max_span = (iteration % 2) ? std::min(n, longspan(rnd)) : n;
			options.max_loop = ((iteration / 2) % 2) ? (longspan(rnd) % 30) : 30;
			options.allow_fft = ((iteration / 4) % 2) == 0;
			options.sequence = MakeRandomRnaSeq(n, rnd());
			options.reference_structure1 = std::string("");
			options.S1 = RandomStructurePK(options.sequence, rnd());
			options.max_dim1 = ComputeMaxHammingDistance(options.sequence, options.S1, options.max_span, options.max_loop);
			const auto ans1 = SimpleMcCaskill<Floating>(options.sequence, options.param_file_name, options.temperature, options.max_span, options.max_loop);
			const auto ans2 = ComputeRintD1Dim<Comp>(options);

			double max_scale = 0.0;
			for (const auto c : ans2) max_scale = std::max<double>(max_scale, double(c.log_scale));

			double mori_sum = 0.0;
			for (const auto c : ans2) {
				if (c.IsZero())continue;
				assert(abs(c.real - 1.0) * exp(c.log_scale - max_scale) < 0.0000000001);
				assert(abs(c.imag) * exp(c.log_scale - max_scale) < 0.0000000001);
				mori_sum += exp(double(c.log_scale));
			}
			assert(abs(mori_sum - double(ans1.second)) < std::max<double>(mori_sum, double(ans1.second)) * 0.0000000001);
		}
	}
}

}

