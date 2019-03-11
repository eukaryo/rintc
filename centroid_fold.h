/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#ifndef RINTDWR_CENTROID_FOLD_H_
#define RINTDWR_CENTROID_FOLD_H_

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

namespace rintdwr {

std::string GetCentroidFoldHamadaBook(
	const std::string& sequence,
	const std::vector<std::vector<Floating>>& bpp_matrix,
	const Floating& gamma,
	const int max_loop);

std::string GetCentroidFoldMcCaskill(
	const std::string& sequence,
	const int max_span,
	const std::vector<std::vector<Floating>>& bpp_matrix,
	const Floating& gamma,
	const int max_loop);

std::vector<std::string> GetCentroidFoldForEachHammingDistance(
	const std::string& sequence,
	const std::vector<std::vector<int>>& S,
	const int max_dim,
	const int max_span,
	const std::vector<std::vector<Floating>>& bpp_matrix,
	const Floating& gamma,
	const int max_loop,
	const int distance = -1);

Floating ComputeGain(
	const std::string& structure,
	const std::vector<std::vector<Floating>>& bpp_matrix,
	const Floating& gamma,
	const int max_loop);

std::pair<std::string, Floating> BruteForceGetCentroidFoldMcCaskill(
	const std::string& sequence,
	const int max_span,
	const std::vector<std::vector<Floating>>& bpp_matrix,
	const Floating& gamma,
	const int max_loop);

std::vector<std::pair<std::string, Floating>> BruteForceGetCentroidFoldForEachHammingDistance(
	const std::string& sequence,
	const std::string& reference_structure,
	const int max_dim,
	const int max_span,
	const std::vector<std::vector<Floating>>& bpp_matrix,
	const Floating& gamma,
	const int max_loop);

}

#endif//RINTDWR_CENTROID_FOLD_H_

