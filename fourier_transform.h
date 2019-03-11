/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#ifndef RINTDWR_FOURIER_TRANSFORM_H_
#define RINTDWR_FOURIER_TRANSFORM_H_

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

std::vector<WideComplexNumber<IntervalVar>>FourierTransform(const std::vector<WideComplexNumber<IntervalVar>>& s, const bool allow_fft);
std::vector<WideComplexNumber<Floating>>FourierTransform(const std::vector<WideComplexNumber<Floating>>& s, const bool allow_fft);
std::vector<UsualComplexNumber<IntervalVar>>FourierTransform(const std::vector<UsualComplexNumber<IntervalVar>>& s, const bool allow_fft);
std::vector<UsualComplexNumber<Floating>>FourierTransform(const std::vector<UsualComplexNumber<Floating>>& s, const bool allow_fft);

}

#endif//RINTDWR_FOURIER_TRANSFORM_H_