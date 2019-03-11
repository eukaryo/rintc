/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#ifndef RINTDWR_TEST_H_
#define RINTDWR_TEST_H_

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

void TestMaxHamming(const int num);
void TestSampling(const int num);
void TestSimpleMcCaskillWide(const int num);
void TestRintD1Dim(const int num);
void TestRintD2Dim(const int num);
void TestRintW1Dim(const int num);
void TestRintW2Dim(const int num);
void TestCentroidFold(const int num);
void TestHagioNonFourier(const int num);
void TestMaxHammingPK(const int num);
void TestRintD1DimPK(const int num);

}

#endif//RINTDWR_TEST_H_