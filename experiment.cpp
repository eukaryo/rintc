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

//#ifdef _MSC_VER
//// 実験的な実装を使う
//#include <experimental/filesystem>
//namespace fs = std::experimental::filesystem;
//#else
//// ファイルシステムをサポートしている
//#include <filesystem>
//namespace fs = std::filesystem;
//#endif

//#include <experimental/filesystem>
//namespace fs = std::experimental::filesystem;


#include"complex_number.h"
#include"misc.h"
#include"rintx.h"
#include"mccaskill_1990.h"
#include"centroid_fold.h"
#include"real_logsumexp.h"

#include"sample_mccaskill.h"
#include"parameter.h"


namespace rintdwr {

static std::vector<std::string> Split(const std::string& str, const char delimiter) {
	std::vector<std::string> ans;
	std::stringstream s(str);
	std::string tmp;
	while (std::getline(s, tmp, delimiter))ans.push_back(tmp);
	return ans;
}

static std::vector<std::vector<std::string>>ReadS151RfamDataset() {
	std::vector<std::vector<std::string>>dataset;
	for (int i = 0; i < 510; ++i) {
		for (const std::string j : std::vector<std::string>{ "A","B" }) {
			const std::string filename =
				"S151RfamDataset/RF00" +
				std::to_string(i / 100) +
				std::to_string((i / 10) % 10) +
				std::to_string(i % 10) +
				"_" + j + ".bpseq";
			std::ifstream ifs(filename);
			if (ifs.fail())continue;
			std::vector<std::string>data{ filename,"","" };
			while (1) {
				std::string line;
				if (!std::getline(ifs, line))break;
				if (line[line.size() - 1] == '\r') {
					line.erase(line.size() - 1);
				}
				const auto v = Split(line, ' ');
				const int pos = stoi(v[0]);
				data[1] += v[1];
				const int bp = stoi(v[2]);
				if (bp == 0)data[2] += ".";
				else if (bp < pos)data[2] += ")";
				else if (bp > pos)data[2] += "(";
				else assert(0);
			}
			dataset.push_back(data);
		}
	}
	assert(int(dataset.size()) == 151);

	//int n = 0;
	//for (int i = 0; i < dataset.size(); i++) {
	//	std::cout << i << " " << dataset[i][1].size() << std::endl;
	//	if (dataset[n][1].size() > dataset[i][1].size())n = i;
	//}
	//std::cout << "shortest length = " << dataset[n][1].size() << std::endl;
	return dataset;
}

static std::string ReadFasta1(const std::string& filename) {

	std::ifstream ifs(filename);
	if (ifs.fail())return std::string("");
	std::string data("");
	while (1) {
		std::string line;
		if (!std::getline(ifs, line))break;
		if (line[0] == '>')continue;
		if (line[line.size() - 1] == '\r') {
			line.erase(line.size() - 1);
		}
		data += line;
	}

	for (int i = 0; i < data.size(); ++i) {
		switch (data[i]) {
		case 'a':data[i] = 'A'; break;
		case 'c':data[i] = 'C'; break;
		case 'g':data[i] = 'G'; break;
		case 't':data[i] = 'U'; break;
		case 'T':data[i] = 'U'; break;
		case 'u':data[i] = 'U'; break;
		}
	}

	return data;

}

static std::pair<std::string,std::string> Read_tRNA() {
	//gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi19/genes/tRNA-Ala-AGC-1-1.html
	std::string seq("GGGGGUAUAGCUCAGUGGUAGAGCGCGUGCUUAGCAUGCACGAGGUCCUGGGUUCGAUCCCCAGUACCUCCA");
	std::string str("(((((((..((((.......)))).(((((.......))))).....(((((.......)))))))))))).");
	return make_pair(seq, str);
}

static std::vector<std::vector<int>>ReadStructure(const std::string& filename, const std::string& sequence) {

	const int n = sequence.size();
	const std::string basepair("AU CG GC GU UA UG");
	std::vector<std::vector<int>>ans(n + 1, std::vector<int>(n + 1, 0));
	std::vector<int>paired(n, 0);

	std::ifstream ifs(filename);
	if (ifs.fail())assert(0);// return ans;
	int count = 0;
	while (1) {
		std::string line;
		if (!std::getline(ifs, line))break;
		if (line[line.size() - 1] == '\r') {
			line.erase(line.size() - 1);
		}
		if (line[0] == '>')continue;
		const auto str = Split(line, ' ');
		if (str.size() != 2)continue;
		const int bp0 = stoi(str[0]);
		const int bp1 = stoi(str[1]);
		assert(bp0 + TURN + 1 <= bp1);
		std::string b("XX");
		b[0] = sequence[bp0];
		b[1] = sequence[bp1];
		assert(basepair.find(b) != std::string::npos);
		assert(paired[bp0]++ == 0);
		assert(paired[bp1]++ == 0);
		ans[bp0 + 1][bp1 + 1]++;
		count++;
	}
	std::cout << count << std::endl;

	return ans;
}


int ComputeBestMaxSpan(
	const std::string& sequence,
	const std::string param_file_name,
	const double temperature,
	const int max_loop,
	const int sample_amount,
	const double threshold_parameter) {

	const int n = int(sequence.size());

	//const auto samples = SampleMcCaskillEnergyAware(sequence, "Turner2004", 37.0, sample_amount, n, max_loop, false, 123456).first;
	//std::vector<std::vector<int>>bp_count(n + 1, std::vector<int>(n + 1, 0));
	//std::vector<int>dist_count(n + 1, 0);
	//for (const std::string structure : samples) {
	//	const auto s = VerifyAndParseStructure(structure, sequence, n, max_loop);
	//	std::vector<int>dist_count_tmp(n + 1, 0);
	//	for (int i = 1; i <= n; ++i)for (int j = i + 1; j <= n; ++j)if(s[i][j]){
	//		bp_count[i][j]++;
	//		dist_count_tmp[j - i] = 1;
	//	}
	//	for (int i = 1; i <= n; ++i)dist_count[i] += dist_count_tmp[i];
	//}
	//std::vector<int>meta_dist_count(n + 1, 0);
	//for (int i = 1; i <= n; ++i)for (int j = i + 1; j <= n; ++j) {
	//	meta_dist_count[j - i] = std::max(meta_dist_count[j - i], bp_count[i][j]);
	//}

	//for (int i = n - 20; i <= n; ++i) {
	//	std::cout << i << " " << meta_dist_count[i] << " " << dist_count[i] << std::endl;
	//}

	//for (int i = n; i >= 0; --i) {
	//	if (double(dist_count[i]) / double(sample_amount) >= 1.0 - threshold_parameter) {
	//		return i;
	//	}
	//}

	return n;

}


void ComputationTimeExperiment1() {
	//random seq

	typedef WideComplexNumber<Floating> Comp;

	RintX1DOptions options;
	options.temperature = 37.0;
	options.param_file_name = std::string("Turner2004");
	options.max_loop = 30;
	options.allow_fft = false;
	const Floating gamma = Floating(1.0);

	std::vector<std::string>sequences;
	for (int i = 100; i <= 800; i += 100)for (int j = 0; j < 1; ++j) {
		sequences.push_back(MakeRandomRnaSeq(i, i + j));
	}

	std::vector<std::vector<int>>result;

	for (int m = 100; m <= 100; m *= 2) {
		for (int i = 0; i < int(sequences.size()); ++i) {
			options.sequence = sequences[i];
			const int length = int(options.sequence.size());
			options.max_span = std::min(length, m);
			options.allow_fft = false;

			std::cout << "start" << std::endl;
			const auto start = std::chrono::system_clock::now();

			const auto bppm = SimpleMcCaskill<Floating>(options.sequence, options.param_file_name, options.temperature, options.max_span, options.max_loop).first;
			options.reference_structure1 = GetCentroidFoldMcCaskill(options.sequence, options.max_span, bppm, gamma, options.max_loop);
			options.S1 = VerifyAndParseStructure(options.reference_structure1, options.sequence, options.max_span, options.max_loop);
			options.max_dim1 = ComputeMaxHammingDistance(options.sequence, options.S1, options.max_span, options.max_loop);
			const auto ans2 = ComputeRintW1Dim<Comp>(options);

			const auto end = std::chrono::system_clock::now();
			std::cout << "end " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;

			const int time = int(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
			std::cout << "result: random seq, i = " << i
				<< ", length = " << length
				<< ", max_span = " << options.max_span
				<< ", time = " << time
				<< " ms" << std::endl;

			result.push_back(std::vector<int>{length, options.max_span, time});
		}
	}
	std::ofstream outputfile("result_random_seq_computation_time.txt");
	for (const auto r : result) {
		outputfile << r[0] << " " << r[1] << " " << r[2] << std::endl;
	}
	outputfile.close();
}
void ComputationTimeExperiment2() {
	//S151Rfam Dataset, maxspan=100

	typedef WideComplexNumber<Floating> Comp;

	RintX1DOptions options;
	options.temperature = 37.0;
	options.param_file_name = std::string("Turner2004");
	options.max_loop = 30;
	options.allow_fft = false;
	const Floating gamma = Floating(1.0);

	const auto dataset = ReadS151RfamDataset();

	std::vector<std::string>sequences;
	for (const auto d : dataset) {
		sequences.push_back(d[1]);
	}

	std::vector<std::vector<int>>result;

	for (int m = 100; m <= 100; m *= 2) {
		for (int i = 0; i < int(sequences.size()); ++i) {
			options.sequence = sequences[i];
			const int length = int(options.sequence.size());
			options.max_span = std::min(length, m);
			options.allow_fft = false;

			std::cout << "start" << std::endl;
			const auto start = std::chrono::system_clock::now();

			const auto bppm = SimpleMcCaskill<Floating>(options.sequence, options.param_file_name, options.temperature, options.max_span, options.max_loop).first;
			options.reference_structure1 = GetCentroidFoldMcCaskill(options.sequence, options.max_span, bppm, gamma, options.max_loop);
			options.S1 = VerifyAndParseStructure(options.reference_structure1, options.sequence, options.max_span, options.max_loop);
			options.max_dim1 = ComputeMaxHammingDistance(options.sequence, options.S1, options.max_span, options.max_loop);
			const auto ans2 = ComputeRintW1Dim<Comp>(options);

			const auto end = std::chrono::system_clock::now();
			std::cout << "end " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;

			const int time = int(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
			std::cout << "result: S151Rfam, i = " << i
				<< ", length = " << length
				<< ", max_span = " << options.max_span
				<< ", time = " << time
				<< " ms" << std::endl;

			result.push_back(std::vector<int>{length, options.max_span, time});
		}
	}
	std::ofstream outputfile("result_S151Rfam_seq_computation_time100.txt");
	for (const auto r : result) {
		outputfile << r[0] << " " << r[1] << " " << r[2] << std::endl;
	}
	outputfile.close();
}
void ComputationTimeExperiment3() {
	//S151Rfam Dataset, maxspan=10000 (i.e. without maxspan constraint)

	typedef WideComplexNumber<Floating> Comp;

	RintX1DOptions options;
	options.temperature = 37.0;
	options.param_file_name = std::string("Turner2004");
	options.max_loop = 30;
	options.allow_fft = false;
	const Floating gamma = Floating(1.0);

	const auto dataset = ReadS151RfamDataset();

	std::vector<std::string>sequences;
	for (const auto d : dataset) {
		sequences.push_back(d[1]);
	}

	std::vector<std::vector<int>>result;

	for (int m = 10000; m <= 10000; m *= 2) {
		for (int i = 0; i < int(sequences.size()); ++i) {
			options.sequence = sequences[i];
			const int length = int(options.sequence.size());
			options.max_span = std::min(length, m);
			options.allow_fft = false;

			std::cout << "start" << std::endl;
			const auto start = std::chrono::system_clock::now();

			const auto bppm = SimpleMcCaskill<Floating>(options.sequence, options.param_file_name, options.temperature, options.max_span, options.max_loop).first;
			options.reference_structure1 = GetCentroidFoldMcCaskill(options.sequence, options.max_span, bppm, gamma, options.max_loop);
			options.S1 = VerifyAndParseStructure(options.reference_structure1, options.sequence, options.max_span, options.max_loop);
			options.max_dim1 = ComputeMaxHammingDistance(options.sequence, options.S1, options.max_span, options.max_loop);
			const auto ans2 = ComputeRintW1Dim<Comp>(options);

			const auto end = std::chrono::system_clock::now();
			std::cout << "end " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;

			const int time = int(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
			std::cout << "result: S151Rfam, i = " << i
				<< ", length = " << length
				<< ", max_span = " << options.max_span
				<< ", time = " << time
				<< " ms" << std::endl;

			result.push_back(std::vector<int>{length, options.max_span, time});
		}
	}
	std::ofstream outputfile("result_S151Rfam_seq_computation_time10000.txt");
	for (const auto r : result) {
		outputfile << r[0] << " " << r[1] << " " << r[2] << std::endl;
	}
	outputfile.close();
}
void ComputationTimeExperiment4() {
	//S151Rfam Dataset, maxspan=10000 (i.e. without maxspan constraint)

	typedef WideComplexNumber<Floating> Comp;

	RintX1DOptions options;
	options.temperature = 37.0;
	options.param_file_name = std::string("Turner2004");
	options.max_loop = 30;
	options.allow_fft = false;
	const Floating gamma = Floating(1.0);

	const auto dataset = ReadS151RfamDataset();

	std::vector<std::string>sequences;
	for (const auto d : dataset) {
		sequences.push_back(d[1]);
	}

	std::vector<std::vector<int>>result;

	for (int m = 10000; m <= 10000; m *= 2) {
		for (int i = 0; i < int(sequences.size()); ++i) {

			options.sequence = sequences[i];
			const int length = int(options.sequence.size());
			options.max_span = std::min(length, m);
			options.allow_fft = false;

			if (length > 100)continue;

			std::cout << "start" << std::endl;
			const auto start = std::chrono::system_clock::now();

			const auto bppm = SimpleMcCaskill<Floating>(options.sequence, options.param_file_name, options.temperature, options.max_span, options.max_loop).first;
			options.reference_structure1 = GetCentroidFoldMcCaskill(options.sequence, options.max_span, bppm, gamma, options.max_loop);
			options.S1 = VerifyAndParseStructure(options.reference_structure1, options.sequence, options.max_span, options.max_loop);
			options.max_dim1 = ComputeMaxHammingDistance(options.sequence, options.S1, options.max_span, options.max_loop);
			const auto ans2 = ComputeRintW1DimNonFourier<Comp>(options);

			const auto end = std::chrono::system_clock::now();
			std::cout << "end " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;

			const int time = int(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
			std::cout << "result: S151Rfam, i = " << i
				<< ", length = " << length
				<< ", max_span = " << options.max_span
				<< ", time = " << time
				<< " ms" << std::endl;

			result.push_back(std::vector<int>{length, options.max_span, time});
		}
	}
	std::ofstream outputfile("result_S151Rfam_seq_computation_timeNF10000.txt");
	for (const auto r : result) {
		outputfile << r[0] << " " << r[1] << " " << r[2] << std::endl;
	}
	outputfile.close();
}

void SuppMaxDistExperiment1() {
	//max hamming distance, S151Rfam Dataset

	const Floating gamma = Floating(1.0);
	const auto dataset = ReadS151RfamDataset();

	std::vector<std::pair<int, int>>result;
	for (int i = 0; i < dataset.size(); ++i) {
		const auto sequence = dataset[i][1];
		const int length = int(sequence.size());
		const int max_span = length;
		const int max_loop = 30;

		std::cout << i << " 1" << std::endl;
		const auto bppm = SimpleMcCaskillWide(sequence, "Turner2004", 37.0, max_span, max_loop).first;
		std::cout << i << " 2" << std::endl;
		const auto reference_structure1 = GetCentroidFoldMcCaskill(sequence, max_span, bppm, gamma, max_loop);
		std::cout << i << " 3" << std::endl;
		const auto S1 = VerifyAndParseStructure(reference_structure1, sequence, max_span, max_loop);
		std::cout << i << " 4" << std::endl;
		const int max_dim1 = ComputeMaxHammingDistance(sequence, S1, max_span, max_loop);
		std::cout << i << " " << length << " " << max_dim1 << std::endl;

		result.push_back(std::make_pair(length, max_dim1));
	}

	std::ofstream outputfile("result_max_hamming_distance_s151.txt");
	for (const auto r : result) {
		outputfile << r.first << " " << r.second << std::endl;
	}
	outputfile.close();
}
void SuppMaxDistExperiment2() {
	//max hamming distance, random sequence, various length

	const Floating gamma = Floating(1.0);

	std::vector<std::pair<int, int>>result;
	for (int i = 400; i >= 50; --i) {

		const auto sequence = MakeRandomRnaSeq(i, i);
		const int length = int(sequence.size());
		const int max_span = length;
		const int max_loop = 30;

		std::cout << i << " 1" << std::endl;
		const auto bppm = SimpleMcCaskillWide(sequence, "Turner2004", 37.0, max_span, max_loop).first;
		std::cout << i << " 2" << std::endl;
		const auto reference_structure1 = GetCentroidFoldMcCaskill(sequence, max_span, bppm, gamma, max_loop);
		std::cout << i << " 3" << std::endl;
		const auto S1 = VerifyAndParseStructure(reference_structure1, sequence, max_span, max_loop);
		std::cout << i << " 4" << std::endl;
		const int max_dim1 = ComputeMaxHammingDistance(sequence, S1, max_span, max_loop);
		std::cout << i << " " << length << " " << max_dim1 << std::endl;

		result.push_back(std::make_pair(length, max_dim1));
	}

	std::ofstream outputfile("result_max_hamming_distance_random.txt");
	for (const auto r : result) {
		outputfile << r.first << " " << r.second << std::endl;
	}
	outputfile.close();
}
void SuppMaxDistExperiment3() {
	//max hamming distance, random sequence, length=200

	const Floating gamma = Floating(1.0);

	std::vector<std::pair<int, int>>result;
	for (int i = 0; i < 1000; ++i) {

		const auto sequence = MakeRandomRnaSeq(200, 1000+i);
		const int length = int(sequence.size());
		const int max_span = length;
		const int max_loop = 30;

		std::cout << i << " / 1000" << std::endl;
		std::cout << i << " 1" << std::endl;
		const auto bppm = SimpleMcCaskillWide(sequence, "Turner2004", 37.0, max_span, max_loop).first;
		std::cout << i << " 2" << std::endl;
		const auto reference_structure1 = GetCentroidFoldMcCaskill(sequence, max_span, bppm, gamma, max_loop);
		std::cout << i << " 3" << std::endl;
		const auto S1 = VerifyAndParseStructure(reference_structure1, sequence, max_span, max_loop);
		std::cout << i << " 4" << std::endl;
		const int max_dim1 = ComputeMaxHammingDistance(sequence, S1, max_span, max_loop);
		std::cout << i << " " << length << " " << max_dim1 << std::endl;

		result.push_back(std::make_pair(length, max_dim1));
	}

	std::ofstream outputfile("result_max_hamming_distance_random200.txt");
	for (const auto r : result) {
		outputfile << r.first << " " << r.second << std::endl;
	}
	outputfile.close();
}

void AccuracyExperiment1(const int num, const bool perform_non_fourier) {
	//numerical accuracy experiment, S151Rfam Dataset

	typedef WideComplexNumber<IntervalVar> Comp;
	const Floating gamma = Floating(1.0);

	RintX1DOptions options;
	options.temperature = 37.0;
	options.param_file_name = std::string("Turner2004");
	options.max_loop = 30;

	const auto s151 = ReadS151RfamDataset();
	std::cout << num << " 1" << std::endl;
	const int n = s151[num][1].size();
	options.max_span = std::min(100, n);
	const auto bppm = SimpleMcCaskillWide(s151[num][1], "Turner2004", 37.0, options.max_span, 30).first;
	std::cout << num << " 2" << std::endl;
	const auto reference_structure = GetCentroidFoldMcCaskill(s151[num][1], options.max_span, bppm, gamma, 30);
	std::cout << num << " 3" << std::endl;
	const auto data = make_pair(s151[num][1], reference_structure);//Read_tRNA();

	options.sequence = data.first;
	options.reference_structure1 = data.second;
	options.S1 = VerifyAndParseStructure(options.reference_structure1, options.sequence, options.max_span, options.max_loop);
	options.max_dim1 = ComputeMaxHammingDistance(options.sequence, options.S1, options.max_span, options.max_loop);

	if(perform_non_fourier)
	{
		std::cout << "start: non fourier" << std::endl;
		const auto ans1 = ComputeRintW1DimNonFourier<Comp>(options);
		const auto ans2 = RintW1DimToBppm(ans1.first, ans1.second);
		std::cout << "end: non fourier" << std::endl;

		std::vector<std::pair<Floating, Floating>>result;
		for (int d = 0; d <= options.max_dim1; ++d) {
			const Floating x = mid(ans2.first[d]);
			for (int i = 1; i <= n; ++i)for (int j = i + 1; j <= n && j - i <= options.max_span; ++j) {
				const Floating y = ans2.second[d][i][j - i].upper() - ans2.second[d][i][j - i].lower();
				//std::cout << x << " " << y << std::endl;
				result.push_back(std::make_pair(x, y));
			}
		}
		std::ofstream outputfile(std::string("result_accuracy_non_fourier_s151_") + std::to_string(num) + std::string(".txt"));
		for (const auto r : result) {
			outputfile << r.first << " " << r.second << std::endl;
		}
		outputfile.close();
	}

	{
		options.allow_fft = false;

		std::cout << "start: dft" << std::endl;
		const auto ans1 = ComputeRintW1Dim<Comp>(options);
		const auto ans2 = RintW1DimToBppm(ans1.first, ans1.second);
		std::cout << "end: dft" << std::endl;

		std::vector<std::pair<Floating, Floating>>result;
		for (int d = 0; d <= options.max_dim1; ++d) {
			const Floating x = mid(ans2.first[d]);
			for (int i = 1; i <= n; ++i)for (int j = i + 1; j <= n && j - i <= options.max_span; ++j) {
				const Floating y = ans2.second[d][i][j - i].upper() - ans2.second[d][i][j - i].lower();
				//std::cout << x << " " << y << std::endl;
				result.push_back(std::make_pair(x, y));
			}
		}
		std::ofstream outputfile(std::string("result_accuracy_dft_s151_") + std::to_string(num) + std::string(".txt"));
		for (const auto r : result) {
			outputfile << r.first << " " << r.second << std::endl;
		}
		outputfile.close();
	}

	{
		options.allow_fft = true;

		std::cout << "start: fft" << std::endl;
		const auto ans1 = ComputeRintW1Dim<Comp>(options);
		const auto ans2 = RintW1DimToBppm(ans1.first, ans1.second);
		std::cout << "end:fft" << std::endl;

		std::vector<std::pair<Floating, Floating>>result;
		for (int d = 0; d <= options.max_dim1; ++d) {
			const Floating x = mid(ans2.first[d]);
			for (int i = 1; i <= n; ++i)for (int j = i + 1; j <= n && j - i <= options.max_span; ++j) {
				const Floating y = ans2.second[d][i][j - i].upper() - ans2.second[d][i][j - i].lower();
				//std::cout << x << " " << y << std::endl;
				result.push_back(std::make_pair(x, y));
			}
		}
		std::ofstream outputfile(std::string("result_accuracy_fft_s151_") + std::to_string(num) + std::string(".txt"));
		for (const auto r : result) {
			outputfile << r.first << " " << r.second << std::endl;
		}
		outputfile.close();
	}
}
void AccuracyExperiment2() {
	for (int i = 0; i < 151; ++i) {
		AccuracyExperiment1(i, false);
	}
}


void AccuracyGlitchVisualization() {

	typedef WideComplexNumber<Floating> Comp;
	const Floating gamma = Floating(1.0);
	RintX1DOptions options;
	options.temperature = 37.0;
	options.param_file_name = std::string("Turner2004");
	options.max_loop = 30;

	const auto s151 = ReadS151RfamDataset();
	options.sequence = s151[6][1];
	options.max_span = options.sequence.size();

	options.allow_fft = true;

	std::cout << " 1" << std::endl;
	const auto bppm = SimpleMcCaskillWide(options.sequence, "Turner2004", 37.0, options.max_span, options.max_loop).first;
	std::cout << " 2" << std::endl;
	options.reference_structure1 = GetCentroidFoldMcCaskill(options.sequence, options.max_span, bppm, gamma, options.max_loop);
	for (int i = 0; i < options.reference_structure1.size(); ++i)options.reference_structure1[i] = '.';
	std::cout << " 3" << std::endl;
	options.S1 = VerifyAndParseStructure(options.reference_structure1, options.sequence, options.max_span, options.max_loop);
	std::cout << " 4" << std::endl;
	options.max_dim1 = ComputeMaxHammingDistance(options.sequence, options.S1, options.max_span, options.max_loop);
	const auto ans11 = ComputeRintW1Dim<Comp>(options);
	const auto ans21 = ComputeRintW1DimNonFourier<Comp>(options);
	const auto ans12 = RintW1DimToBppm(ans11.first, ans11.second);
	const auto ans22 = RintW1DimToBppm(ans21.first, ans21.second);

	std::vector<std::string> represent1, represent2;
	for (int i = 0; i <= options.max_dim1; ++i) {
		represent1.push_back(GetCentroidFoldForEachHammingDistance(options.sequence, options.S1, options.max_dim1, options.max_span, ans12.second[i], gamma, options.max_loop, i)[0]);
		represent2.push_back(GetCentroidFoldForEachHammingDistance(options.sequence, options.S1, options.max_dim1, options.max_span, ans22.second[i], gamma, options.max_loop, i)[0]);
		//represent1.push_back(GetCentroidFoldMcCaskill(options.sequence, options.max_span, ans12.second[i], gamma, options.max_loop));
		//represent2.push_back(GetCentroidFoldMcCaskill(options.sequence, options.max_span, ans22.second[i], gamma, options.max_loop));
	}

	for (int i = 0; i <= options.max_dim1; ++i) {
		std::cout << i << std::endl << ans22.first[i] << std::endl << represent1[i] << std::endl << represent2[i] << std::endl << std::endl;
	}

	return;


}

void HeatResistanceExperiment(const std::string filename, double threshold) {
	//CentroidFoldでリファレンス取ってHagio→取り直してMori

	//temp <- [37, 55]

	const std::string rrna = ReadFasta1(filename);
	const Floating gamma = Floating(1.0);
	const int span_default = 100;

	RintX1DOptions options;

	std::cout << "start: filename = " << filename << std::endl;
	std::cout << "prepare: credibility limit" << std::endl;

	options.param_file_name = std::string("Turner2004");
	options.max_loop = 30;
	options.allow_fft = false;
	options.temperature = 37.0;

	options.sequence = rrna;
	options.max_span = std::min(int(options.sequence.size()), span_default);
	std::cout << "1" << std::endl;
	const auto bppm = SimpleMcCaskillWide(options.sequence, options.param_file_name, options.temperature, options.max_span, options.max_loop).first;
	std::cout << "2" << std::endl;
	const std::string old_ref = GetCentroidFoldMcCaskill(options.sequence, options.max_span, bppm, gamma, options.max_loop);
	options.reference_structure1 = old_ref;
	std::cout << "3" << std::endl;
	options.S1 = VerifyAndParseStructure(options.reference_structure1, options.sequence, options.max_span, options.max_loop);
	std::cout << "4" << std::endl;
	options.max_dim1 = ComputeMaxHammingDistance(options.sequence, options.S1, options.max_span, options.max_loop);

	//{
	//	std::vector<std::vector<int>>result;
	//	std::cout << "start: credibility limit(old ref)" << std::endl;
	//	for (int i = 37; i <= 55; i += 2) {
	//		options.temperature = double(i);

	//		const int n = int(options.sequence.size());
	//		const auto tmp = ComputeRintD1Dim<WideComplexNumber<Floating>>(options);
	//		const auto mori_plot = RegularizeRintD1Dim(tmp);
	//		Floating sum = Floating(0.0);
	//		int CL = -1;
	//		for (int j = 0; j <= options.max_dim1; ++j) {
	//			sum += mori_plot[j];
	//			if (sum >= threshold) {
	//				CL = j;
	//				break;
	//			}
	//		}

	//		std::cout << i << " " << CL << " (old ref)" << std::endl;
	//		result.push_back(std::vector<int>{i, CL});

	//		std::ofstream outputfile(std::string("result_mori_plot(old ref)_temp_") + std::to_string(i) + std::string("_") + filename + std::string(".txt"));
	//		for (int j = 0; j <= options.max_dim1; ++j) {
	//			outputfile << j << " " << mori_plot[j] << std::endl;
	//		}
	//		outputfile.close();

	//	}
	//	std::cout << "end: credibility limit(old ref)" << std::endl;

	//	std::ofstream outputfile(std::string("result_credibility_limit(old ref)_") + filename + std::string(".txt"));
	//	outputfile << old_ref << std::endl;
	//	for (const auto r : result) {
	//		outputfile << r[0] << " " << r[1] << std::endl;
	//	}
	//	outputfile.close();
	//}

	std::cout << "start: Hagiotool for making better reference" << std::endl;
	const auto hagio = ComputeRintW1Dim<WideComplexNumber<Floating>>(options);
	std::cout << "ToBppm..." << std::endl;
	const auto hagio_result = RintW1DimToBppm(hagio.first, hagio.second);
	std::cout << "end: Hagiotool for making better reference" << std::endl;
	//山のBPPMを求める

	//vector<float>を入力して、連続する部分列で総和が0.1以上になるような最短のものの左端のインデックスを出力する関数を作る。
	auto func1 = [&](const double sum) {
		int ansL = -1, ansLen = hagio_result.first.size();
		for (int i = 0; i < hagio_result.first.size(); ++i) {
			double s = 0.0;
			for (int j = i; j < hagio_result.first.size(); ++j) {
				s += hagio_result.first[j];
				if (s >= sum) {
					if (ansLen > j - i + 1) {
						ansLen = j - i + 1;
						ansL = i;
					}
					break;
				}
			}
		}
		assert(ansL != -1);
		return ansL;
	};

	//左端のインデックスを出力する。
	auto func2 = [&](const double sum, const int ansL) {
		double actual_sum = 0.0;
		int ansR = -1;
		for (int i = ansL; i < hagio_result.first.size(); ++i) {
			actual_sum += hagio_result.first[i];
			if (actual_sum >= sum) {
				ansR = i;
				break;
			}
		}
		return ansR;
	};

	//新しいリファレンスを作るためのBPPMを求める。
	auto func3 = [&](const int ansL, const int ansR) {

		double actual_sum = 0.0;
		for (int i = ansL; i <= ansR; ++i)actual_sum += hagio_result.first[i];
		auto ans_bppm = hagio_result.second[ansL];
		for (int a = 0; a < ans_bppm.size(); ++a)for (int b = 0; b < ans_bppm[a].size(); ++b)ans_bppm[a][b] = 0.0;

		for (int i = ansL; i <= ansR; ++i) {
			for (int a = 0; a < ans_bppm.size(); ++a) {
				for (int b = 0; b < ans_bppm[a].size(); ++b) {
					ans_bppm[a][b] += hagio_result.second[ansL][a][b] * (hagio_result.first[ansL] / actual_sum);
				}
			}
		}
		return ans_bppm;
	};

	const int ansL = func1(0.1);
	const int ansR = func2(0.1, ansL);
	const auto new_bppm = func3(ansL, ansR);
	std::cout << "4" << std::endl;
	const std::string new_ref = GetCentroidFoldMcCaskill(options.sequence, options.max_span, new_bppm, gamma, options.max_loop);

	//int maxindex = 0;
	//for (int i = 1; i <= options.max_dim1; ++i)if (hagio_result.first[maxindex] < hagio_result.first[i])maxindex = i;
	//std::cout << maxindex << std::endl;
	//const auto peak_bppm = hagio_result.second[maxindex];
	//const std::string new_ref = GetCentroidFoldForEachHammingDistance(options.sequence, options.S1, options.max_dim1, options.max_span, peak_bppm, gamma, options.max_loop, maxindex)[0];
	//assert(ComputeHammingDistance(old_ref, new_ref) == maxindex);

	options.reference_structure1 = new_ref;
	std::cout << "5" << std::endl;
	options.S1 = VerifyAndParseStructure(options.reference_structure1, options.sequence, options.max_span, options.max_loop);
	std::cout << "6" << std::endl;
	options.max_dim1 = ComputeMaxHammingDistance(options.sequence, options.S1, options.max_span, options.max_loop);

	{
		std::vector<std::vector<int>>result;
		std::cout << "start: credibility limit(new ref)" << std::endl;
		for (int i = 37; i <= 55; i += 2) {
			options.temperature = double(i);

			const int n = int(options.sequence.size());
			const auto tmp = ComputeRintD1Dim<WideComplexNumber<Floating>>(options);
			const auto mori_plot = RegularizeRintD1Dim(tmp);
			Floating sum = Floating(0.0);
			int CL = -1;
			for (int j = 0; j <= options.max_dim1; ++j) {
				sum += mori_plot[j];
				if (sum >= threshold) {
					CL = j;
					break;
				}
			}

			std::cout << i << " " << CL << " (new ref)" << std::endl;
			result.push_back(std::vector<int>{i, CL});

			std::ofstream outputfile(std::string("result_mori_plot(new ref)_temp_") + std::to_string(i) + std::string("_") + filename + std::string(".txt"));
			for (int j = 0; j <= options.max_dim1; ++j) {
				outputfile << j << " " << mori_plot[j] << std::endl;
			}
			outputfile.close();

		}
		std::cout << "end: credibility limit(new ref)" << std::endl;

		std::ofstream outputfile(std::string("result_credibility_limit(new ref)_") + filename + std::string(".txt"));
		outputfile << ansL << std::endl;
		outputfile << ansR << std::endl;
		outputfile << new_ref << std::endl;
		for (const auto r : result) {
			outputfile << r[0] << " " << r[1] << std::endl;
		}
		outputfile.close();
	}
}

void HeatResistanceExperimentPK(const std::string sequencefilename, const std::string structurefilename, double threshold) {

	//temp <- [37, 55]

	const std::string rrna = ReadFasta1(sequencefilename);
	const Floating gamma = Floating(1.0);
	const int span_default = 100;

	RintX1DOptions options;

	std::cout << "start: filename = " << sequencefilename << std::endl;
	std::cout << "prepare: credibility limit" << std::endl;

	options.param_file_name = std::string("Turner2004");
	options.max_loop = 30;
	options.allow_fft = false;
	options.temperature = 37.0;

	options.sequence = rrna;
	options.max_span = std::min(int(options.sequence.size()), span_default);
	std::cout << "1" << std::endl;
	std::cout << "2" << std::endl;
	options.reference_structure1 = std::string("");
	std::cout << "3" << std::endl;
	options.S1 = ReadStructure(structurefilename, options.sequence);
	std::cout << "4" << std::endl;
	options.max_dim1 = ComputeMaxHammingDistance(options.sequence, options.S1, options.max_span, options.max_loop);

	std::vector<std::vector<int>>result;
	std::cout << "start: credibility limit" << std::endl;
	for (int i = 37; i <= 55; i += 2) {
		options.temperature = double(i);

		const int n = int(options.sequence.size());
		const auto tmp = ComputeRintD1Dim<WideComplexNumber<Floating>>(options);
		const auto mori_plot = RegularizeRintD1Dim(tmp);
		Floating sum = Floating(0.0);
		int CL = -1;
		for (int j = 0; j <= options.max_dim1; ++j) {
			sum += mori_plot[j];
			if (sum >= threshold) {
				CL = j;
				break;
			}
		}

		std::cout << i << " " << CL << " " << std::endl;
		result.push_back(std::vector<int>{i, CL});

		std::ofstream outputfile(std::string("result_refPK_mori_plot_temp_") + std::to_string(i) + std::string("_") + sequencefilename + std::string(".txt"));
		for (int j = 0; j <= options.max_dim1; ++j) {
			outputfile << j << " " << mori_plot[j] << std::endl;
		}
		outputfile.close();

	}
	std::cout << "end: credibility limit" << std::endl;

	std::ofstream outputfile(std::string("result_refPK_credibility_limit_") + sequencefilename + std::string(".txt"));
	for (const auto r : result) {
		outputfile << r[0] << " " << r[1] << std::endl;
	}
	outputfile.close();

}

void RrnaAccuracyExperiment(const std::string filename) {
	//CentroidFoldでリファレンス取ってHagi

	//temp 37

	const std::string rrna = ReadFasta1(filename);
	const Floating gamma = Floating(1.0);
	const int span_default = 100;

	RintX1DOptions options;

	std::cout << "start: filename = " << filename << std::endl;
	std::cout << "prepare: credibility limit" << std::endl;

	options.param_file_name = std::string("Turner2004");
	options.max_loop = 30;
	options.allow_fft = false;
	options.temperature = 37.0;

	options.sequence = rrna;
	options.max_span = std::min(int(options.sequence.size()), span_default);
	std::cout << "1" << std::endl;
	const auto bppm = SimpleMcCaskillWide(options.sequence, options.param_file_name, options.temperature, options.max_span, options.max_loop).first;
	std::cout << "2" << std::endl;
	const std::string old_ref = GetCentroidFoldMcCaskill(options.sequence, options.max_span, bppm, gamma, options.max_loop);
	options.reference_structure1 = old_ref;
	std::cout << "3" << std::endl;
	options.S1 = VerifyAndParseStructure(options.reference_structure1, options.sequence, options.max_span, options.max_loop);
	std::cout << "4" << std::endl;
	options.max_dim1 = ComputeMaxHammingDistance(options.sequence, options.S1, options.max_span, options.max_loop);

	typedef WideComplexNumber<IntervalVar> Comp;
	const int n = options.sequence.size();

	{
		options.allow_fft = false;

		std::cout << "start: dft" << std::endl;
		const auto ans1 = ComputeRintW1Dim<Comp>(options);
		const auto ans2 = RintW1DimToBppm(ans1.first, ans1.second);
		std::cout << "end: dft" << std::endl;

		std::vector<std::pair<Floating, Floating>>result;
		for (int d = 0; d <= options.max_dim1; ++d) {
			const Floating x = mid(ans2.first[d]);
			for (int i = 1; i <= n; ++i)for (int j = i + 1; j <= n && j - i <= options.max_span; ++j) {
				const Floating y = ans2.second[d][i][j - i].upper() - ans2.second[d][i][j - i].lower();
				//std::cout << x << " " << y << std::endl;
				result.push_back(std::make_pair(x, y));
			}
		}
		std::ofstream outputfile(std::string("result_accuracy_dft_rrna_") + filename + std::string(".txt"));
		for (const auto r : result) {
			outputfile << r.first << " " << r.second << std::endl;
		}
		outputfile.close();
	}

	{
		options.allow_fft = true;

		std::cout << "start: fft" << std::endl;
		const auto ans1 = ComputeRintW1Dim<Comp>(options);
		const auto ans2 = RintW1DimToBppm(ans1.first, ans1.second);
		std::cout << "end:fft" << std::endl;

		std::vector<std::pair<Floating, Floating>>result;
		for (int d = 0; d <= options.max_dim1; ++d) {
			const Floating x = mid(ans2.first[d]);
			for (int i = 1; i <= n; ++i)for (int j = i + 1; j <= n && j - i <= options.max_span; ++j) {
				const Floating y = ans2.second[d][i][j - i].upper() - ans2.second[d][i][j - i].lower();
				//std::cout << x << " " << y << std::endl;
				result.push_back(std::make_pair(x, y));
			}
		}
		std::ofstream outputfile(std::string("result_accuracy_fft_rrna_") + filename + std::string(".txt"));
		for (const auto r : result) {
			outputfile << r.first << " " << r.second << std::endl;
		}
		outputfile.close();
	}

}


}
