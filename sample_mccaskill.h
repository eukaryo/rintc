/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#ifndef RINTDWR_SAMPLE_MCCASKILL_H_
#define RINTDWR_SAMPLE_MCCASKILL_H_

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

//#include <boost/multiprecision/cpp_int.hpp>

#include"real_logsumexp.h"

namespace rintdwr {

//typedef boost::multiprecision::cpp_int Bigint;
typedef int64_t Bigint;
typedef WideRealNumber<double> WideFloating;

std::pair<std::vector<std::string>, Bigint>SampleMcCaskillUniform(
	const std::string sequence,
	const int sample_amount,
	const int max_span,
	const int max_loop,
	const bool do_debug,
	const int seed);

std::pair<std::vector<std::string>, WideFloating>SampleMcCaskillEnergyAware(
	const std::string sequence,
	const std::string param_file_name,
	const double temperature,
	const int sample_amount,
	const int max_span,
	const int max_loop,
	const bool do_debug,
	const int seed);
}


#endif//RINTDWR_SAMPLE_MCCASKILL_H_