/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#ifndef RINTDWR_MCCASKILL_1990_H_
#define RINTDWR_MCCASKILL_1990_H_

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

#include"interval_type.h"
#include"real_logsumexp.h"

namespace rintdwr {

template<typename RealScalar>std::pair<std::vector<std::vector<RealScalar>>, RealScalar>SimpleMcCaskill(
	const std::string sequence,
	const std::string param_file_name,
	const double temperature,
	const int max_span,
	const int max_loop);

void OutputBppm(
	const std::string sequence,
	const std::string param_file_name,
	const double temperature,
	const int max_span,
	const int max_loop);

std::pair<std::vector<std::vector<Floating>>, Floating>SimpleMcCaskillWide(
	const std::string sequence,
	const std::string param_file_name,
	const double temperature,
	const int max_span,
	const int max_loop);
}


#endif//RINTDWR_MCCASKILL_1990_H_