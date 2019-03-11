/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#ifndef RINTDWR_RINTX_H_
#define RINTDWR_RINTX_H_

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

#include"complex_number.h"
#include"interval_type.h"

namespace rintdwr {

//Google C++ Style Guide:
//> For functions that have several configuration options,
//> consider defining a single class or struct to hold all the options,
//> and pass an instance of that.

//コードを大体書いた後でこの記述に気付いた。後付け的に↓のstructを作った。
//(次に取り組むプロジェクトでは最初から気を付けて設計しようと思うが、今はさておく)

struct RintX1DOptions {
	std::string sequence;
	std::string reference_structure1;
	std::vector<std::vector<int>> S1;//VerifyAndParseStructure(reference_structure1, sequence, max_span, max_loop)
	int max_dim1;
	std::string param_file_name;
	double temperature;
	int max_span;
	int max_loop;
	bool allow_fft;

	void Dump()const {
		std::cout << "sequence  : " << sequence << std::endl;
		std::cout << "structure1: " << reference_structure1 << std::endl;
		std::cout << "max_dim1  : " << max_dim1 << std::endl;
		std::cout << "param_file: " << param_file_name << std::endl;
		std::cout << "temp      : " << temperature << std::endl;
		std::cout << "max_span  : " << max_span << std::endl;
		std::cout << "max_loop  : " << max_loop << std::endl;
	}
};

struct RintX2DOptions : RintX1DOptions {
	std::string reference_structure2;
	std::vector<std::vector<int>> S2;//VerifyAndParseStructure(reference_structure2, sequence, max_span, max_loop)
	int max_dim2;

	void Dump()const {
		std::cout << "sequence  : " << sequence << std::endl;
		std::cout << "structure1: " << reference_structure1 << std::endl;
		std::cout << "structure1: " << reference_structure2 << std::endl;
		std::cout << "max_dim1  : " << max_dim1 << std::endl;
		std::cout << "max_dim2  : " << max_dim2 << std::endl;
		std::cout << "param_file: " << param_file_name << std::endl;
		std::cout << "temp      : " << temperature << std::endl;
		std::cout << "max_span  : " << max_span << std::endl;
		std::cout << "max_loop  : " << max_loop << std::endl;
	}

};

template<typename Comp>std::vector<Comp>ComputeRintD1Dim(const RintX1DOptions& options);
template<typename Comp>std::vector<std::vector<Comp>>ComputeRintD2Dim(const RintX2DOptions& options);
template<typename Comp>std::pair<std::vector<Comp>, std::vector<std::vector<std::vector<Comp>>>> ComputeRintW1Dim(const RintX1DOptions& options);
template<typename Comp>std::pair<std::vector<Comp>, std::vector<std::vector<std::vector<Comp>>>> ComputeRintW1DimNonFourier(const RintX1DOptions& options);
template<typename Comp>std::pair<std::vector<std::vector<Comp>>, std::vector<std::vector<std::vector<std::vector<Comp>>>>> ComputeRintW2Dim(const RintX2DOptions& options);

std::vector<double> BruteForceRintD1Dim(const RintX1DOptions& options);
std::vector<double> BruteForceRintD1DimPK(const RintX1DOptions& options);
std::vector<std::vector<double>> BruteForceRintD2Dim(const RintX2DOptions& options);
std::pair<std::vector<double>, std::vector<std::vector<std::vector<double>>>> BruteForceRintW1Dim(const RintX1DOptions& options);
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<std::vector<std::vector<double>>>>> BruteForceRintW2Dim(const RintX2DOptions& options);

std::pair<std::vector<Floating>, std::vector<std::vector<std::vector<Floating>>>>RintW1DimToBppm(
	const std::vector<WideComplexNumber<Floating>> z,
	const std::vector<std::vector<std::vector<WideComplexNumber<Floating>>>> p);

std::pair<std::vector<IntervalVar>, std::vector<std::vector<std::vector<IntervalVar>>>>RintW1DimToBppm(
	const std::vector<WideComplexNumber<IntervalVar>> z,
	const std::vector<std::vector<std::vector<WideComplexNumber<IntervalVar>>>> p);

std::vector<Floating>RegularizeRintD1Dim(const std::vector<WideComplexNumber<Floating>>& result);
std::vector<std::vector<Floating>>RegularizeRintD2Dim(const std::vector<std::vector<WideComplexNumber<Floating>>>& result);

std::pair<std::vector<std::vector<Floating>>, std::vector<std::vector<std::vector<std::vector<Floating>>>>>RintW2DimToBppm(
	const std::vector<std::vector<WideComplexNumber<Floating>>> z,
	const std::vector<std::vector<std::vector<std::vector<WideComplexNumber<Floating>>>>> p);

int ComputeCredibilityLimit(const RintX1DOptions& options, const Floating threshold);

std::vector<std::string> ComputeRintW1Dim_1(
	const RintX1DOptions& options,
	const int dim);

std::vector<std::string> ComputeRintW1Dim_2(
	const RintX1DOptions& options,
	const std::vector<WideComplexNumber<IntervalVar>>& zeta,
	const std::vector<std::vector<WideComplexNumber<IntervalVar>>>& prob,
	const int num,
	const int mod);

}

#endif//RINTDWR_RINTX_H_