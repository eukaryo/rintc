/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#include"real_logsumexp.h"

namespace rintdwr {

std::vector<std::vector<Floating>> ConvertBppmToUsualReal(const std::vector<std::vector<WideRealNumber<Floating>>>& bppm) {
	
	const int x = int(bppm.size());
	const int y = int(bppm[0].size());

	std::vector<std::vector<Floating>> ans(x, std::vector<Floating>(y, Floating(0.0)));

	for (int i = 0; i < x; ++i)for (int j = 0; j < y; ++j) {
		ans[i][j] = bppm[i][j].ToUsualReal();
	}

	return ans;
}

}

