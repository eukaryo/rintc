/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#include"mccaskill_1990.h"

#include"interval_type.h"
#include"real_logsumexp.h"
#include"parameter.h"
#include"misc.h"

namespace rintdwr {

//explicit instantiation

template std::pair<std::vector<std::vector<Floating>>, Floating>SimpleMcCaskill(
	const std::string sequence,
	const std::string param_file_name,
	const double temperature,
	const int max_span,
	const int max_loop);
template std::pair<std::vector<std::vector<IntervalVar>>, IntervalVar>SimpleMcCaskill(
	const std::string sequence,
	const std::string param_file_name,
	const double temperature,
	const int max_span,
	const int max_loop);

//function definition

template<typename RealScalar>std::pair<std::vector<std::vector<RealScalar>>, RealScalar>SimpleMcCaskill(
	const std::string sequence,
	const std::string param_file_name,
	const double temperature,
	const int max_span,
	const int max_loop) {
	//Compute [McCaskill, 1990] and return (base_pairing_probability_matrix, inside_partition_functuin).

	parasor_param::InitializeParameter(param_file_name, temperature);

	const int n = int(sequence.size());

	typedef std::pair<bool, RealScalar>MemReal;

	std::vector<MemReal>Z((n + 1) * 2, std::make_pair(false, RealScalar(0.0)));
	std::vector<MemReal>Z1((n + 1) * (max_span + 1), std::make_pair(false, RealScalar(0.0)));
	std::vector<MemReal>Zb((n + 1) * (max_span + 1), std::make_pair(false, RealScalar(0.0)));
	std::vector<MemReal>Zm((n + 1) * (max_span + 1), std::make_pair(false, RealScalar(0.0)));
	std::vector<MemReal>Zm1((n + 1) * (max_span + 1), std::make_pair(false, RealScalar(0.0)));
	std::vector<MemReal>Wb((n + 1) * (max_span + 1), std::make_pair(false, RealScalar(0.0)));

	const auto at = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n && j - i <= max_span);
		return i * (max_span + 1) + (j - i);
	};

	std::function<RealScalar(int, int)> GetZ;
	std::function<RealScalar(int, int)> GetZ1;
	std::function<RealScalar(int, int)> GetZb;
	std::function<RealScalar(int, int)> GetZm;
	std::function<RealScalar(int, int)> GetZm1;
	std::function<RealScalar(int, int)> GetWb;

	GetZ = [&](const int i, const int j) {
		assert(1 <= i && (i <= j || i == j + 1) && j <= n);
		if (i == j)return RealScalar(1.0);
		if (i == j + 1)return RealScalar(1.0);
		assert(i == 1 || j == n);
		const int index = (i == 1) ? j : (i + n + 1);
		assert(1 <= index && index < (n + 1) * 2);
		if (Z[index].first)return Z[index].second;

		//RealScalar ans = RealScalar(1.0);
		//for (int k = i; k <= j; ++k)ans += GetZ(i, k - 1) * GetZ1(k, j);
		//return (Z[i][j] = std::make_pair(true, ans)).second;

		RealScalar ans = RealScalar(0.0);
		if (i == 1) {

			//no base pair
			ans += RealScalar(1.0);

			//there is only one outermost base pair (i, *)
			ans += GetZ1(i, j);

			//there are one or more outermost base pairs. The rightmost outermost base pair is (k + 1, *)
			for (int k = i; k <= j - 1; ++k) {
				ans += GetZ(i, k) * GetZ1(k + 1, j);
			}
		}
		else {

			//i is not paired
			ans += GetZ(i + 1, j);

			//(i,k) form a base pair
			for (int k = i + TURN + 1; k <= j && k - i <= max_span; ++k) {
				const int type = parasor_param::GetPairType(sequence[i - 1], sequence[k - 1]);
				if (type == 0)continue;
				ans += GetZb(i, k)
					* exp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (k - 1) + 1, true, sequence)))
					* GetZ(k + 1, j);
			}
		}
		return (Z[index] = std::make_pair(true, ans)).second;

	};
	GetZ1 = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n);
		if (i == j)return RealScalar(0.0);

		if ((j - i) > max_span) {
			return GetZ1(i, i + max_span);
		}

		if (Z1[at(i, j)].first)return Z1[at(i, j)].second;

		RealScalar ans = RealScalar(0.0);
		const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
		if (type != 0) {
			ans += GetZb(i, j)
				* exp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, true, sequence)));
		}
		ans += GetZ1(i, j - 1);
		return (Z1[at(i, j)] = std::make_pair(true, ans)).second;
	};
	GetZb = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n);
		if (i == j)return RealScalar(0.0);

		const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
		if (type == 0 || i + TURN >= j || j - i > max_span) {
			return RealScalar(0.0);
		}

		if (Zb[at(i, j)].first)return Zb[at(i, j)].second;

		RealScalar ans = RealScalar(0.0);

		//hairpin loop
		ans += exp(RealScalar(parasor_param::ParHairpinEnergy(i - 1, j - 1, sequence)));

		//internal loop, stem, bulge
		for (int k = i + 1; k <= std::min(i + max_loop + 1, j - TURN - 2); ++k) {
			const int unpaired_base1 = k - i - 1;
			for (int l = std::max(k + TURN + 1, j - 1 - max_loop + unpaired_base1); l < j; ++l) {
				const int type = parasor_param::GetPairType(sequence[k - 1], sequence[l - 1]);
				if (type == 0)continue;
				ans += GetZb(k, l) * exp(RealScalar(parasor_param::ParLoopEnergy(i - 1, j - 1, k - 1, l - 1, sequence)));
			}
		}

		//multi loop
		for (int k = i + TURN + 3; k + TURN + 1 <= j - 1; ++k) {
			ans += GetZm(i + 1, k - 1) * GetZm1(k, j - 1)
				* exp(RealScalar(parasor_param::ParDangling(parasor_param::GetPairTypeReverse(sequence[i - 1], sequence[j - 1]), (j - 1) - 1, (i - 1) + 1, false, sequence)))
				* exp(RealScalar(parasor_param::ParMultiloopClosing()));
		}

		return (Zb[at(i, j)] = std::make_pair(true, ans)).second;
	};
	GetZm = [&](const int i, const int j) {
		assert(1 <= i && (i <= j || i == j + 1) && j <= n);
		if (i == j)return RealScalar(0.0);
		if (i == j + 1)return RealScalar(0.0);
		if (Zm[at(i, j)].first)return Zm[at(i, j)].second;

		RealScalar ans = RealScalar(0.0);
		for (int k = i; k/* + TURN + 1*/ <= j; ++k) {
			ans += (1.0 + GetZm(i, k - 1)) * GetZm1(k, j);
		}
		return (Zm[at(i, j)] = std::make_pair(true, ans)).second;
	};
	GetZm1 = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n);
		if (i == j)return RealScalar(0.0);
		if (Zm1[at(i, j)].first)return Zm1[at(i, j)].second;

		RealScalar ans = RealScalar(0.0);
		const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
		if (type != 0) {
			ans += GetZb(i, j)
				* exp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, false, sequence)))
				* exp(RealScalar(parasor_param::ParMultiloopInternal()));
		}
		ans += GetZm1(i, j - 1);
		return (Zm1[at(i, j)] = std::make_pair(true, ans)).second;
	};
	GetWb = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n);

		const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
		if (type == 0 || i + TURN >= j || j - i > max_span) {
			return RealScalar(0.0);
		}

		if (Wb[at(i, j)].first)return Wb[at(i, j)].second;

		RealScalar ans = RealScalar(0.0);

		ans += GetZ(1, i - 1) * GetZ(j + 1, n) * exp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, true, sequence)));

//		for (int h = 1; h < i; ++h) {
		for (int h = std::max(1, i - max_span + TURN + 2); h < i; ++h) {
//			for (int l = j + 1; l <= n; ++l) {
			for (int l = j + 1; l <= n && (l - h) <= max_span; ++l) {
						const int rtype = parasor_param::GetPairTypeReverse(sequence[h - 1], sequence[l - 1]);
				if (rtype == 0)continue;

				const int u = (i - h - 1) + (l - j - 1);
				if (u <= max_loop) {
					//internal loop, bulge, stem
					ans += GetWb(h, l) * exp(RealScalar(parasor_param::ParLoopEnergy(h - 1, l - 1, i - 1, j - 1, sequence)));
				}

				//(h,l)‚Å•Â‚¶‚ç‚ê‚émultiloop
				ans += exp(RealScalar(parasor_param::ParMultiloopClosing()))
					* exp(RealScalar(parasor_param::ParMultiloopInternal()))
					* GetWb(h, l)
					* exp(RealScalar(parasor_param::ParDangling(rtype, (l - 1) - 1, (h - 1) + 1, false, sequence)))
					* exp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, false, sequence)))
					* (
						GetZm(h + 1, i - 1)
						+ GetZm(j + 1, l - 1)
						+ GetZm(h + 1, i - 1) * GetZm(j + 1, l - 1)
						);
			}
		}
		return (Wb[at(i, j)] = std::make_pair(true, ans)).second;
	};

	auto bppm = std::vector<std::vector<RealScalar>>(n + 1, std::vector<RealScalar>(max_span + 1, RealScalar(0.0)));
	for (int i = 1; i <= n; ++i)for (int j = i; j <= n && j - i <= max_span; ++j) {
		bppm[i][j - i] = GetZb(i, j) * GetWb(i, j) / GetZ(1, n);
	}

	return std::make_pair(bppm, GetZ(1, n));
}

std::pair<std::vector<std::vector<Floating>>, Floating>SimpleMcCaskillWide(
	const std::string sequence,
	const std::string param_file_name,
	const double temperature,
	const int max_span,
	const int max_loop) {
	//Compute [McCaskill, 1990] and return (base_pairing_probability_matrix, log(inside_partition_function)).

	parasor_param::InitializeParameter(param_file_name, temperature);

	const int n = int(sequence.size());

	typedef WideRealNumber<Floating> RealScalar;
	typedef std::pair<bool, WideRealNumber<Floating>>MemReal;

	std::vector<MemReal>Z((n + 1) * 2, std::make_pair(false, RealScalar(0.0)));
	std::vector<MemReal>Z1((n + 1) * (max_span + 1), std::make_pair(false, RealScalar(0.0)));
	std::vector<MemReal>Zb((n + 1) * (max_span + 1), std::make_pair(false, RealScalar(0.0)));
	std::vector<MemReal>Zm((n + 1) * (max_span + 1), std::make_pair(false, RealScalar(0.0)));
	std::vector<MemReal>Zm1((n + 1) * (max_span + 1), std::make_pair(false, RealScalar(0.0)));
	std::vector<MemReal>Wb((n + 1) * (max_span + 1), std::make_pair(false, RealScalar(0.0)));

	const auto at = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n && j - i <= max_span);
		return i * (max_span + 1) + (j - i);
	};

	std::function<RealScalar(int, int)> GetZ;
	std::function<RealScalar(int, int)> GetZ1;
	std::function<RealScalar(int, int)> GetZb;
	std::function<RealScalar(int, int)> GetZm;
	std::function<RealScalar(int, int)> GetZm1;
	std::function<RealScalar(int, int)> GetWb;

	GetZ = [&](const int i, const int j) {
		assert(1 <= i && (i <= j || i == j + 1) && j <= n);
		if (i == j)return RealScalar(1.0);
		if (i == j + 1)return RealScalar(1.0);
		assert(i == 1 || j == n);
		const int index = (i == 1) ? j : (i + n + 1);
		assert(1 <= index && index < (n + 1) * 2);
		if (Z[index].first)return Z[index].second;

		//RealScalar ans = RealScalar(1.0);
		//for (int k = i; k <= j; ++k)ans += GetZ(i, k - 1) * GetZ1(k, j);
		//return (Z[i][j] = std::make_pair(true, ans)).second;

		RealScalar ans = RealScalar(0.0);
		if (i == 1) {

			//no base pair
			ans += RealScalar(1.0);

			//there is only one outermost base pair (i, *)
			ans += GetZ1(i, j);

			//there are one or more outermost base pairs. The rightmost outermost base pair is (k + 1, *)
			for (int k = i; k <= j - 1; ++k) {
				ans += GetZ(i, k) * GetZ1(k + 1, j);
			}
		}
		else {

			//i is not paired
			ans += GetZ(i + 1, j);

			//(i,k) form a base pair
			for (int k = i + TURN + 1; k <= j && k - i <= max_span; ++k) {
				const int type = parasor_param::GetPairType(sequence[i - 1], sequence[k - 1]);
				if (type == 0)continue;
				ans += GetZb(i, k)
					* RealScalar::LogRealToWide(Floating(parasor_param::ParDangling(type, (i - 1) - 1, (k - 1) + 1, true, sequence)))
					* GetZ(k + 1, j);
			}
		}
		return (Z[index] = std::make_pair(true, ans)).second;

	};
	GetZ1 = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n);
		if (i == j)return RealScalar(0.0);

		if ((j - i) > max_span) {
			return GetZ1(i, i + max_span);
		}

		if (Z1[at(i, j)].first)return Z1[at(i, j)].second;

		RealScalar ans = RealScalar(0.0);
		const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
		if (type != 0) {
			ans += GetZb(i, j)
				* RealScalar::LogRealToWide(Floating(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, true, sequence)));
		}
		ans += GetZ1(i, j - 1);
		return (Z1[at(i, j)] = std::make_pair(true, ans)).second;
	};
	GetZb = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n);
		if (i == j)return RealScalar(0.0);

		const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
		if (type == 0 || i + TURN >= j || j - i > max_span) {
			return RealScalar(0.0);
		}

		if (Zb[at(i, j)].first)return Zb[at(i, j)].second;

		RealScalar ans = RealScalar(0.0);

		//hairpin loop
		ans += RealScalar::LogRealToWide(Floating(parasor_param::ParHairpinEnergy(i - 1, j - 1, sequence)));

		//internal loop, stem, bulge
		for (int k = i + 1; k <= std::min(i + max_loop + 1, j - TURN - 2); ++k) {
			const int unpaired_base1 = k - i - 1;
			for (int l = std::max(k + TURN + 1, j - 1 - max_loop + unpaired_base1); l < j; ++l) {
				const int type = parasor_param::GetPairType(sequence[k - 1], sequence[l - 1]);
				if (type == 0)continue;
				ans += GetZb(k, l) * RealScalar::LogRealToWide(Floating(parasor_param::ParLoopEnergy(i - 1, j - 1, k - 1, l - 1, sequence)));
			}
		}

		//multi loop
		for (int k = i + TURN + 3; k + TURN + 1 <= j - 1; ++k) {
			ans += GetZm(i + 1, k - 1) * GetZm1(k, j - 1)
				* RealScalar::LogRealToWide(Floating(parasor_param::ParDangling(parasor_param::GetPairTypeReverse(sequence[i - 1], sequence[j - 1]), (j - 1) - 1, (i - 1) + 1, false, sequence)))
				* RealScalar::LogRealToWide(Floating(parasor_param::ParMultiloopClosing()));
		}

		return (Zb[at(i, j)] = std::make_pair(true, ans)).second;
	};
	GetZm = [&](const int i, const int j) {
		assert(1 <= i && (i <= j || i == j + 1) && j <= n);
		if (i == j)return RealScalar(0.0);
		if (i == j + 1)return RealScalar(0.0);
		if (Zm[at(i, j)].first)return Zm[at(i, j)].second;

		RealScalar ans = RealScalar(0.0);
		for (int k = i; k/* + TURN + 1*/ <= j; ++k) {
			ans += (RealScalar(1.0) + GetZm(i, k - 1)) * GetZm1(k, j);
		}
		return (Zm[at(i, j)] = std::make_pair(true, ans)).second;
	};
	GetZm1 = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n);
		if (i == j)return RealScalar(0.0);
		if (Zm1[at(i, j)].first)return Zm1[at(i, j)].second;

		RealScalar ans = RealScalar(0.0);
		const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
		if (type != 0) {
			ans += GetZb(i, j)
				* RealScalar::LogRealToWide(Floating(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, false, sequence)))
				* RealScalar::LogRealToWide(Floating(parasor_param::ParMultiloopInternal()));
		}
		ans += GetZm1(i, j - 1);
		return (Zm1[at(i, j)] = std::make_pair(true, ans)).second;
	};
	GetWb = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n);

		const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
		if (type == 0 || i + TURN >= j || j - i > max_span) {
			return RealScalar(0.0);
		}

		if (Wb[at(i, j)].first)return Wb[at(i, j)].second;

		RealScalar ans = RealScalar(0.0);

		ans += GetZ(1, i - 1) * GetZ(j + 1, n) * RealScalar::LogRealToWide(Floating(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, true, sequence)));

		//		for (int h = 1; h < i; ++h) {
		for (int h = std::max(1, i - max_span + TURN + 2); h < i; ++h) {
			//			for (int l = j + 1; l <= n; ++l) {
			for (int l = j + 1; l <= n && (l - h) <= max_span; ++l) {
				const int rtype = parasor_param::GetPairTypeReverse(sequence[h - 1], sequence[l - 1]);
				if (rtype == 0)continue;

				const int u = (i - h - 1) + (l - j - 1);
				if (u <= max_loop) {
					//internal loop, bulge, stem
					ans += GetWb(h, l) * RealScalar::LogRealToWide(Floating(parasor_param::ParLoopEnergy(h - 1, l - 1, i - 1, j - 1, sequence)));
				}

				//(h,l)‚Å•Â‚¶‚ç‚ê‚émultiloop
				ans += RealScalar::LogRealToWide(Floating(parasor_param::ParMultiloopClosing()))
					* RealScalar::LogRealToWide(Floating(parasor_param::ParMultiloopInternal()))
					* GetWb(h, l)
					* RealScalar::LogRealToWide(Floating(parasor_param::ParDangling(rtype, (l - 1) - 1, (h - 1) + 1, false, sequence)))
					* RealScalar::LogRealToWide(Floating(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, false, sequence)))
					* (
						GetZm(h + 1, i - 1)
						+ GetZm(j + 1, l - 1)
						+ GetZm(h + 1, i - 1) * GetZm(j + 1, l - 1)
						);
			}
		}
		return (Wb[at(i, j)] = std::make_pair(true, ans)).second;
	};

	auto bppm = std::vector<std::vector<Floating>>(n + 1, std::vector<Floating>(max_span + 1, Floating(0.0)));
	for (int i = 1; i <= n; ++i)for (int j = i; j <= n && j - i <= max_span; ++j) {
		bppm[i][j - i] = (GetZb(i, j) * GetWb(i, j) / GetZ(1, n)).ToUsualReal();
	}

	return std::make_pair(bppm, GetZ(1, n).log_scale);
}


void OutputBppm(
	const std::string sequence,
	const std::string param_file_name,
	const double temperature,
	const int max_span,
	const int max_loop) {

	//show BPPM

	const auto bppm = SimpleMcCaskill<Floating>(sequence, param_file_name, temperature, max_span, max_loop).first;

	for (int i = 1; i < int(sequence.size()); ++i) {
		for (int j = 1; j <= int(sequence.size()) && j - i <= max_span; ++j) {
			if (j != 1)std::cout << ",";
			std::cout << bppm[i][j - i];
		}
		std::cout << std::endl;
	}
}

}
