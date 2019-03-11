/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#ifndef RINTDWR_MISC_H_
#define RINTDWR_MISC_H_

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

namespace rintdwr {

constexpr int TURN = 3;

int Ceiling2Power(int h);
std::string MakeRandomRnaSeq(const int length, const int seed);
int ComputeMaxLoop(const std::string& structure);
std::vector<std::vector<int>> VerifyAndParseStructure(const std::string& structure, const std::string& sequence, const int max_span, const int max_loop);
std::vector<std::vector<int>> ComputePredistanceMatrix(const std::vector<std::vector<int>>& S);
int ComputeMaxHammingDistance(const std::string& sequence, const std::vector<std::vector<int>>& S, const int max_span, const int max_loop);
std::vector<std::string>EnumerateStructures(const std::string& sequence, const int max_span, const int max_loop);
std::vector<std::vector<int>>RandomStructurePK(const std::string& sequence, const int seed);
int ComputeHammingDistance(const std::string& structure1, const std::string& structure2);
int ComputeHammingDistance(const std::vector<std::vector<int>>& structure1, const std::vector<std::vector<int>>& structure2);
std::string ComputeStructuralContext(const std::string& structure, const int pos);
double EvalSpecificStructure(const std::string& sequence, const std::string& structure);
double EvalSpecificStructure(const std::string& sequence, const std::vector<std::vector<int>>& structure);
std::string MatrixToDotNotation(const std::vector<std::vector<int>>& structure);


}

#endif//RINTDWR_MISC_H_