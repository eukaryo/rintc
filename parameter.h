/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#ifndef RINTDWR_PARAMETER_H_
#define RINTDWR_PARAMETER_H_

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

namespace parasor_param {

int GetBaseType(const char c);
int GetPairType(const char c1, const char c2);
int GetPairTypeReverse(const char c1, const char c2);

double ParMultiloopClosing();
double ParMultiloopInternal();
double ParDangling(const int type, const int five, const int three, const bool ext, const std::string& sequence);
double ParHairpinEnergy(const int i, const int j, const std::string& sequence);
double ParLoopEnergy(const int i, const int j, const int p, const int q, const std::string& sequence);
void InitializeParameter(const std::string& parameter_file_name, const double temperature);

extern bool counting;

}

}


#endif//RINTDWR_PARAMETER_H_