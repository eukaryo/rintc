/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#include"sample_mccaskill.h"

#include<random>

#include <boost/random.hpp>

#include"parameter.h"
#include"misc.h"

namespace rintdwr {

std::pair<std::vector<std::string>, Bigint>SampleMcCaskillUniform(
	const std::string sequence,
	const int sample_amount,
	const int max_span,
	const int max_loop,
	const bool do_debug,
	const int seed) {

	//単なるMcCaskill型DPをやって、内側分配関数のstochastic backtrackをしてサンプリングする。
	//エネルギーモデルは導入しない。全ての二次構造から一様にサンプリングする。
	//サンプリングしたRNA構造たちと、内側分配関数の値を返す。

	const int n = int(sequence.size());

	typedef std::pair<bool, Bigint>MemBigint;

	std::vector<MemBigint>Z((n + 1) * 2, std::make_pair(false, Bigint(0)));
	std::vector<MemBigint>Z1((n + 1) * (max_span + 1), std::make_pair(false, Bigint(0)));
	std::vector<MemBigint>Zb((n + 1) * (max_span + 1), std::make_pair(false, Bigint(0)));
	std::vector<MemBigint>Zm((n + 1) * (max_span + 1), std::make_pair(false, Bigint(0)));

	const auto at = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n && j - i <= max_span);
		return i * (max_span + 1) + (j - i);
	};

	std::function<Bigint(int, int)> GetZ;
	std::function<Bigint(int, int)> GetZ1;
	std::function<Bigint(int, int)> GetZb;
	std::function<Bigint(int, int)> GetZm;

	GetZ = [&](const int i, const int j) {
		assert(1 <= i && (i <= j || i == j + 1) && j <= n);
		if (i == j)return Bigint(1);
		if (i == j + 1)return Bigint(1);
		assert(i == 1 || j == n);
		const int index = (i == 1) ? j : (i + n + 1);
		assert(1 <= index && index < (n + 1) * 2);
		if (Z[index].first)return Z[index].second;

		Bigint ans = Bigint(0);
		//for (int k = i; k <= j - 1; ++k)ans += GetZ(i, k - 1) * GetZ1(k, j);
		//return (Z[i][j] = std::make_pair(true, ans)).second;

		if (i == 1) {

			//no base pair
			ans += Bigint(1);

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
				ans += GetZb(i, k) * GetZ(k + 1, j);
			}
		}

		return (Z[index] = std::make_pair(true, ans)).second;
	};
	GetZ1 = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n);
		if (i == j)return Bigint(0);

		if ((j - i) > max_span) {
			return GetZ1(i, i + max_span);
		}
		if (Z1[at(i, j)].first)return Z1[at(i, j)].second;

		Bigint ans = GetZb(i, j) + GetZ1(i, j - 1);
		return (Z1[at(i, j)] = std::make_pair(true, ans)).second;
	};
	GetZb = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n);
		if (i == j)return Bigint(0);

		const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
		if (type == 0 || i + TURN >= j || j - i > max_span) {
			return Bigint(0);
		}

		if (Zb[at(i, j)].first)return Zb[at(i, j)].second;

		Bigint ans = Bigint(0);

		//hairpin loop
		ans += Bigint(1);

		//internal loop, stem, bulge
		for (int k = i + 1; k <= std::min(i + max_loop + 1, j - TURN - 2); ++k) {
			const int unpaired_base1 = k - i - 1;
			for (int l = std::max(k + TURN + 1, j - 1 - max_loop + unpaired_base1); l < j; ++l) {
				const int type = parasor_param::GetPairType(sequence[k - 1], sequence[l - 1]);
				if (type == 0)continue;
				ans += GetZb(k, l);
			}
		}

		//multi loop
		for (int k = i + TURN + 3; k + TURN + 1 <= j - 1; ++k) {
			ans += GetZm(i + 1, k - 1) * GetZ1(k, j - 1);
		}

		return (Zb[at(i, j)] = std::make_pair(true, ans)).second;
	};
	GetZm = [&](const int i, const int j) {
		assert(1 <= i && (i <= j || i == j + 1) && j <= n);
		if (i == j)return Bigint(0);
		if (i == j + 1)return Bigint(0);
		if (Zm[at(i, j)].first)return Zm[at(i, j)].second;

		Bigint ans = Bigint(0);
		for (int k = i; k + TURN + 1 <= j; ++k) {
			ans += (Bigint(1) + GetZm(i, k - 1)) * GetZ1(k, j);
		}
		return (Zm[at(i, j)] = std::make_pair(true, ans)).second;
	};

	const Bigint total = GetZ(1, n);
	//boost::random::mt19937_64 rnd(seed);
	//boost::random::uniform_int_distribution<Bigint>random_structure(Bigint(0), total - Bigint(1));
	std::mt19937_64 rnd(seed);
	std::uniform_int_distribution<Bigint>random_structure(Bigint(0), total - Bigint(1));

	const auto GetRandomStructure = [&]() {

		std::function<std::string(int, int, Bigint)> TrackZ;
		std::function<std::string(int, int, Bigint)> TrackZ1;
		std::function<std::string(int, int, Bigint)> TrackZb;
		std::function<std::string(int, int, Bigint)> TrackZm;

		const auto EmptyStructure = [](const int i, const int j) {
			std::string ans("");
			for (int k = i; k <= j; k++)ans += std::string(".");
			return ans;
		};

		TrackZ = [&](const int i, const int j, Bigint index) {
			assert(1 <= i && (i <= j || i == j + 1) && j <= n);
			if (i == j)return std::string(".");
			if (i == j + 1)return std::string("");
			assert(Bigint(0) <= index && index < GetZ(i, j));

			for (int k = i; k <= j - 1; ++k) {
				const Bigint x = GetZ(i, k - 1);
				const Bigint y = GetZ1(k, j);
				const Bigint weight = x * y;
				if (index < weight)return TrackZ(i, k - 1, index % x) + TrackZ1(k, j, index / x);
				else index -= weight;
			}

			assert(index == Bigint(0));
			return EmptyStructure(i, j);
		};
		TrackZ1 = [&](const int i, const int j, Bigint index) {
			assert(1 <= i && i <= j && j <= n);
			if (i == j)return std::string(".");
			assert(Bigint(0) <= index && index < GetZ1(i, j));

			const Bigint weight = GetZb(i, j);
			if (index < weight)return TrackZb(i, j, index);
			return TrackZ1(i, j - 1, index - weight) + std::string(".");
		};
		TrackZb = [&](const int i, const int j, Bigint index) {
			assert(1 <= i && i <= j && j <= n);
			if (i == j)return std::string(".");
			assert(Bigint(0) <= index && index < GetZb(i, j));

			const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
			if (type == 0 || i + TURN >= j || j - i > max_span) {
				assert(0);
				return EmptyStructure(i, j);
			}

			//internal loop, stem, bulge
			for (int k = i + 1; k <= std::min(i + max_loop + 1, j - TURN - 2); ++k) {
				const int unpaired_base1 = k - i - 1;
				for (int l = std::max(k + TURN + 1, j - 1 - max_loop + unpaired_base1); l < j; ++l) {
					const int type = parasor_param::GetPairType(sequence[k - 1], sequence[l - 1]);
					if (type == 0)continue;
					const Bigint weight_internal = GetZb(k, l);
					if (index < weight_internal) {
						return std::string("(") + EmptyStructure(i + 1, k - 1) + TrackZb(k, l, index) + EmptyStructure(l + 1, j - 1) + std::string(")");
					}
					else index -= weight_internal;
				}
			}

			//multi loop
			for (int k = i + TURN + 3; k + TURN + 1 <= j - 1; ++k) {
				const Bigint x = GetZm(i + 1, k - 1);
				const Bigint y = GetZ1(k, j - 1);
				const Bigint weight_multi = x * y;
				if (index < weight_multi)return std::string("(") + TrackZm(i + 1, k - 1, index % x) + TrackZ1(k, j - 1, index / x) + std::string(")");
				else index -= weight_multi;
			}

			//hairpin loop
			assert(index == Bigint(0));
			return std::string("(") + EmptyStructure(i + 1, j - 1) + std::string(")");
		};
		TrackZm = [&](const int i, const int j, Bigint index) {
			assert(1 <= i && (i <= j || i == j + 1) && j <= n);
			if (i == j)return std::string(".");
			if (i == j + 1)return std::string("");
			assert(Bigint(0) <= index && index < GetZm(i, j));

			for (int k = i; k + TURN + 1 <= j; ++k) {
				const Bigint weight1 = GetZ1(k, j);
				if (index < weight1)return EmptyStructure(i, k - 1) + TrackZ1(k, j, index);
				else index -= weight1;

				const Bigint x = GetZm(i, k - 1);
				const Bigint y = GetZ1(k, j);
				const Bigint weight2 = x * y;
				if (index < weight2)return TrackZm(i, k - 1, index % x) + TrackZ1(k, j, index / x);
				else index -= weight2;
			}

			assert(0);
			return EmptyStructure(i, j);
		};

		return TrackZ(1, n, random_structure(rnd));
	};

	if (do_debug) {
		//debug:
		//(1)dump total number of secondary structures (verify using another tool)
		//(2)Compute BPPM with outside_algorithm and sampling. Those BPPMs must be nearly the same.

		std::vector<MemBigint>Wb((n + 1) * (max_span + 1), std::make_pair(false, Bigint(0)));
		std::function<Bigint(int, int)> GetWb;

		GetWb = [&](const int i, const int j) {
			assert(1 <= i && i <= j && j <= n);

			const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
			if (type == 0 || i + TURN >= j) {
				return Bigint(0);
			}

			if (Wb[at(i, j)].first)return Wb[at(i, j)].second;

			Bigint ans = Bigint(0);

			ans += GetZ(1, i - 1) * GetZ(j + 1, n);

			for (int h = std::max(1, i - max_span + TURN + 2); h < i; ++h) {
				for (int l = j + 1; l <= n && (l - h) <= max_span; ++l) {

					const int rtype = parasor_param::GetPairType(sequence[h - 1], sequence[l - 1]);
					if (rtype == 0)continue;

					const int u = (i - h - 1) + (l - j - 1);
					if (u <= max_loop) {
						//internal loop, bulge, stem
						ans += GetWb(h, l);
					}

					//multiloop
					ans += GetWb(h, l) * (
						GetZm(h + 1, i - 1)
						+ GetZm(j + 1, l - 1)
						+ GetZm(h + 1, i - 1) * GetZm(j + 1, l - 1)
						);
				}
			}
			return (Wb[at(i, j)] = std::make_pair(true, ans)).second;
		};

		auto bppm = std::vector<std::vector<Bigint>>(n + 1, std::vector<Bigint>(max_span + 1, Bigint(0)));

		for (int i = 1; i <= n; ++i)for (int j = i; j <= n && j - i <= max_span; ++j) {
			bppm[i][j - i] = GetZb(i, j) * GetWb(i, j);
		}
		std::vector<std::vector<Bigint>>pseudo_bppm(n + 1, std::vector<Bigint>(max_span + 1, Bigint(0)));

		constexpr Bigint size__ = 100000;
		for (int x = 0; x < size__; x++) {
			const std::string structure = GetRandomStructure();
			const auto p = VerifyAndParseStructure(structure, sequence, max_span, max_loop);
			for (int i = 1; i <= n; ++i)for (int j = i + 1; j <= n && j - i <= max_span; ++j)pseudo_bppm[i][j - i] += Bigint(p[i][j]);
		}
		for (int i = 1; i <= n; ++i)for (int j = i + 1; j <= n && j - i <= max_span; ++j) {

			//diff = abs(pseudo_bppm[i][j - i]/size__ - bppm[i][j - i]/total;
			//assert(diff < 0.01);

			const Bigint diff = std::abs<Bigint>(pseudo_bppm[i][j - i] * total - bppm[i][j - i] * size__);
			if (diff * Bigint(100) >= total * Bigint(size__)) {
				std::cout << std::endl;
				std::cout << i << std::endl;
				std::cout << j << std::endl;
				std::cout << max_loop << std::endl;
				std::cout << diff << std::endl;
				std::cout << pseudo_bppm[i][j - i] << std::endl;
				std::cout << total << std::endl;
				std::cout << bppm[i][j - i] << std::endl;
				assert(0);
			}
		}
	}

	std::vector<std::string>ans;

	for (int i = 0; i < sample_amount; i++) {
		const std::string structure = GetRandomStructure();
		VerifyAndParseStructure(structure, sequence, max_span, max_loop);
		ans.push_back(structure);
	}

	return std::make_pair(ans, total);
}

std::pair<std::vector<std::string>, WideFloating>SampleMcCaskillEnergyAware(
	const std::string sequence,
	const std::string param_file_name,
	const double temperature,
	const int sample_amount,
	const int max_span,
	const int max_loop,
	const bool do_debug,
	const int seed) {

	//単なるMcCaskill型DPをやって、内側分配関数のstochastic backtrackをしてサンプリングする。
	//エネルギーモデルを考慮して、ボツルマン因子に比例する確率でサンプリングする。
	//サンプリングしたRNA構造たちと、内側分配関数の値を返す。

	parasor_param::InitializeParameter(param_file_name, temperature);

	const int n = int(sequence.size());

	typedef std::pair<bool, WideFloating>MemWideFloating;

	std::vector<MemWideFloating>Z((n + 1) * 2, std::make_pair(false, WideFloating(0.0)));
	std::vector<MemWideFloating>Z1((n + 1) * (max_span + 1), std::make_pair(false, WideFloating(0.0)));
	std::vector<MemWideFloating>Zb((n + 1) * (max_span + 1), std::make_pair(false, WideFloating(0.0)));
	std::vector<MemWideFloating>Zm((n + 1) * (max_span + 1), std::make_pair(false, WideFloating(0.0)));
	std::vector<MemWideFloating>Zm1((n + 1) * (max_span + 1), std::make_pair(false, WideFloating(0.0)));

	const auto at = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n && j - i <= max_span);
		return i * (max_span + 1) + (j - i);
	};

	std::function<WideFloating(int, int)> GetZ;
	std::function<WideFloating(int, int)> GetZ1;
	std::function<WideFloating(int, int)> GetZb;
	std::function<WideFloating(int, int)> GetZm;
	std::function<WideFloating(int, int)> GetZm1;

	GetZ = [&](const int i, const int j) {
		assert(1 <= i && (i <= j || i == j + 1) && j <= n);
		if (i == j)return WideFloating(1.0);
		if (i == j + 1)return WideFloating(1.0);
		assert(i == 1 || j == n);
		const int index = (i == 1) ? j : (i + n + 1);
		assert(1 <= index && index < (n + 1) * 2);
		if (Z[index].first)return Z[index].second;

		WideFloating ans = WideFloating(0.0);

		if (i == 1) {

			//no base pair
			ans += WideFloating(1.0);

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
					* GetZ(k + 1, j)
					* WideFloating::LogRealToWide(parasor_param::ParDangling(type, (i - 1) - 1, (k - 1) + 1, true, sequence));
			}
		}

		return (Z[index] = std::make_pair(true, ans)).second;
	};
	GetZ1 = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n);
		if (i == j)return WideFloating(0.0);

		if ((j - i) > max_span) {
			return GetZ1(i, i + max_span);
		}
		if (Z1[at(i, j)].first)return Z1[at(i, j)].second;

		WideFloating ans = WideFloating(0.0);
		const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
		if (type != 0) {
			ans += GetZb(i, j)
				* WideFloating::LogRealToWide(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, true, sequence));
		}
		ans += GetZ1(i, j - 1);
		return (Z1[at(i, j)] = std::make_pair(true, ans)).second;
	};
	GetZb = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n);
		if (i == j)return WideFloating(0.0);

		const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
		if (type == 0 || i + TURN >= j || j - i > max_span) {
			return WideFloating(0.0);
		}

		if (Zb[at(i, j)].first)return Zb[at(i, j)].second;

		WideFloating ans(0.0);

		//hairpin loop
		ans += WideFloating::LogRealToWide(parasor_param::ParHairpinEnergy(i - 1, j - 1, sequence));

		//internal loop, stem, bulge
		for (int k = i + 1; k <= std::min(i + max_loop + 1, j - TURN - 2); ++k) {
			const int unpaired_base1 = k - i - 1;
			for (int l = std::max(k + TURN + 1, j - 1 - max_loop + unpaired_base1); l < j; ++l) {
				const int type = parasor_param::GetPairType(sequence[k - 1], sequence[l - 1]);
				if (type == 0)continue;
				ans += GetZb(k, l)
					* WideFloating::LogRealToWide(parasor_param::ParLoopEnergy(i - 1, j - 1, k - 1, l - 1, sequence));
			}
		}

		//multi loop
		const int rtype = parasor_param::GetPairTypeReverse(sequence[i - 1], sequence[j - 1]);
		for (int k = i + TURN + 3; k + TURN + 1 <= j - 1; ++k) {
			ans += GetZm(i + 1, k - 1)
				* GetZm1(k, j - 1)
				* WideFloating::LogRealToWide(parasor_param::ParDangling(rtype, (j - 1) - 1, (i - 1) + 1, false, sequence))
				* WideFloating::LogRealToWide(parasor_param::ParMultiloopClosing());
		}

		return (Zb[at(i, j)] = std::make_pair(true, ans)).second;
	};
	GetZm = [&](const int i, const int j) {
		assert(1 <= i && (i <= j || i == j + 1) && j <= n);
		if (i == j)return WideFloating(0.0);
		if (i == j + 1)return WideFloating(0.0);
		if (Zm[at(i, j)].first)return Zm[at(i, j)].second;

		WideFloating ans = WideFloating(0.0);
		for (int k = i; k + TURN + 1 <= j; ++k) {
			ans += (WideFloating(1.0) + GetZm(i, k - 1)) * GetZm1(k, j);
		}
		return (Zm[at(i, j)] = std::make_pair(true, ans)).second;
	};
	GetZm1 = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n);
		if (i == j)return WideFloating(0.0);
		if (Zm1[at(i, j)].first)return Zm1[at(i, j)].second;

		WideFloating ans = WideFloating(0.0);
		const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
		if (type != 0) {
			ans += GetZb(i, j)
				* WideFloating::LogRealToWide(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, false, sequence))
				* WideFloating::LogRealToWide(parasor_param::ParMultiloopInternal());
		}
		ans += GetZm1(i, j - 1);

		return (Zm1[at(i, j)] = std::make_pair(true, ans)).second;
	};

	const WideFloating total = GetZ(1, n);
	std::mt19937_64 rnd(seed);
	std::uniform_real_distribution<double>random_number(0.0, 1.0);

	const auto GetRandomStructure = [&]() {

		std::function<std::string(int, int)> TrackZ;
		std::function<std::string(int, int)> TrackZ1;
		std::function<std::string(int, int)> TrackZb;
		std::function<std::string(int, int)> TrackZm;
		std::function<std::string(int, int)> TrackZm1;

		const auto EmptyStructure = [](const int i, const int j) {
			std::string ans("");
			for (int k = i; k <= j; k++)ans += std::string(".");
			return ans;
		};

		TrackZ = [&](const int i, const int j) {
			assert(1 <= i && (i <= j || i == j + 1) && j <= n);
			if (i == j)return std::string(".");
			if (i == j + 1)return std::string("");

			WideFloating partition_function = GetZ(i, j);

			for (int k = i; k <= j - 1; ++k) {
				const WideFloating weight = GetZ(i, k - 1) * GetZ1(k, j);
				const double probability = WideFloating::Probability(weight, partition_function);
				if (random_number(rnd) < probability)return TrackZ(i, k - 1) + TrackZ1(k, j);
				else partition_function -= weight;
			}

			return EmptyStructure(i, j);
		};
		TrackZ1 = [&](const int i, const int j) {
			assert(1 <= i && i <= j && j <= n);
			if (i == j)return std::string(".");

			const WideFloating partition_function = GetZ1(i, j);
			const WideFloating weight = GetZ1(i, j - 1);
			const double probability = WideFloating::Probability(weight, partition_function);
			if (random_number(rnd) < probability)return TrackZ1(i, j - 1) + std::string(".");
			return TrackZb(i, j);
		};
		TrackZb = [&](const int i, const int j) {
			assert(1 <= i && i <= j && j <= n);
			if (i == j)return std::string(".");

			const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
			if (type == 0 || i + TURN >= j || j - i > max_span) {
				assert(0);
				return EmptyStructure(i, j);
			}

			WideFloating partition_function = GetZb(i, j);

			//internal loop, stem, bulge
			for (int k = i + 1; k <= std::min(i + max_loop + 1, j - TURN - 2); ++k) {
				const int unpaired_base1 = k - i - 1;
				for (int l = std::max(k + TURN + 1, j - 1 - max_loop + unpaired_base1); l < j; ++l) {
					const int type = parasor_param::GetPairType(sequence[k - 1], sequence[l - 1]);
					if (type == 0)continue;
					const WideFloating weight_internal =
						GetZb(k, l)
						* WideFloating::LogRealToWide(parasor_param::ParLoopEnergy(i - 1, j - 1, k - 1, l - 1, sequence));
					;
					const double probability_internal = WideFloating::Probability(weight_internal, partition_function);
					if (random_number(rnd) < probability_internal) {
						return std::string("(") + EmptyStructure(i + 1, k - 1) + TrackZb(k, l) + EmptyStructure(l + 1, j - 1) + std::string(")");
					}
					else partition_function -= weight_internal;
				}
			}

			//multi loop
			const int rtype = parasor_param::GetPairTypeReverse(sequence[i - 1], sequence[j - 1]);
			for (int k = i + TURN + 3; k + TURN + 1 <= j - 1; ++k) {
				const WideFloating weight_multi =
					GetZm(i + 1, k - 1)
					* GetZm1(k, j - 1)
					* WideFloating::LogRealToWide(parasor_param::ParDangling(rtype, (j - 1) - 1, (i - 1) + 1, false, sequence))
					* WideFloating::LogRealToWide(parasor_param::ParMultiloopClosing());
				const double probability_multi = WideFloating::Probability(weight_multi, partition_function);
				if (random_number(rnd) < probability_multi) {
					return std::string("(") + TrackZm(i + 1, k - 1) + TrackZm1(k, j - 1) + std::string(")");
				}
				else partition_function -= weight_multi;
			}

			//hairpin loop
			return std::string("(") + EmptyStructure(i + 1, j - 1) + std::string(")");
		};
		TrackZm = [&](const int i, const int j) {
			assert(1 <= i && (i <= j || i == j + 1) && j <= n);
			if (i == j)return std::string(".");
			if (i == j + 1)return std::string("");

			WideFloating partition_function = GetZm(i, j);

			for (int k = i; k + TURN + 1 <= j; ++k) {
				const WideFloating weight1 = GetZm1(k, j);
				const double probability1 = WideFloating::Probability(weight1, partition_function);
				if (random_number(rnd) < probability1)return EmptyStructure(i, k - 1) + TrackZm1(k, j);
				else partition_function -= weight1;

				const WideFloating weight2 = GetZm(i, k - 1) * GetZm1(k, j);
				const double probability2 = WideFloating::Probability(weight2, partition_function);
				if (random_number(rnd) < probability2)return TrackZm(i, k - 1) + TrackZm1(k, j);
				else partition_function -= weight2;
			}

			//数値誤差が無ければ、ここに来る可能性は0%である。
			assert(0);
			return EmptyStructure(i, j);
		};
		TrackZm1 = [&](const int i, const int j) {
			assert(1 <= i && i <= j && j <= n);
			if (i == j)return std::string(".");

			const WideFloating partition_function = GetZm1(i, j);
			const WideFloating weight = GetZm1(i, j - 1);
			const double probability = WideFloating::Probability(weight, partition_function);
			if (random_number(rnd) < probability)return TrackZm1(i, j - 1) + std::string(".");
			return TrackZb(i, j);
		};

		return TrackZ(1, n);
	};

	if (do_debug) {
		//debug:
		//(1)dump total number of secondary structures (verify using another tool)
		//(2)Compute BPPM with outside_algorithm and sampling. Those BPPMs must be nearly the same.

		std::vector<MemWideFloating>Wb((n + 1) * (max_span + 1), std::make_pair(false, WideFloating(0.0)));
		std::function<WideFloating(int, int)> GetWb;

		GetWb = [&](const int i, const int j) {
			assert(1 <= i && i <= j && j <= n);

			const int type = parasor_param::GetPairType(sequence[i - 1], sequence[j - 1]);
			if (type == 0 || i + TURN >= j) {
				return WideFloating(0.0);
			}

			if (Wb[at(i, j)].first)return Wb[at(i, j)].second;

			WideFloating ans = WideFloating(0.0);

			ans += GetZ(1, i - 1)
				* GetZ(j + 1, n)
				* WideFloating::LogRealToWide(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, true, sequence));

			for (int h = std::max(1, i - max_span + TURN + 2); h < i; ++h) {
				for (int l = j + 1; l <= n && (l - h) <= max_span; ++l) {

					const int rtype = parasor_param::GetPairType(sequence[h - 1], sequence[l - 1]);
					if (rtype == 0)continue;

					const int u = (i - h - 1) + (l - j - 1);
					if (u <= max_loop) {
						//internal loop, bulge, stem
						ans += GetWb(h, l)
							* WideFloating::LogRealToWide(parasor_param::ParLoopEnergy(h - 1, l - 1, i - 1, j - 1, sequence));
					}

					//multiloop
					ans += WideFloating::LogRealToWide(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, false, sequence))
						* WideFloating::LogRealToWide(parasor_param::ParMultiloopInternal())
						* GetWb(h, l)
						* WideFloating::LogRealToWide(parasor_param::ParDangling(rtype, (l - 1) - 1, (h - 1) + 1, false, sequence))
						* WideFloating::LogRealToWide(parasor_param::ParMultiloopClosing())
						* (
							GetZm(h + 1, i - 1)
							+ GetZm(j + 1, l - 1)
							+ GetZm(h + 1, i - 1) * GetZm(j + 1, l - 1)
							);
				}
			}
			return (Wb[at(i, j)] = std::make_pair(true, ans)).second;
		};

		auto bppm = std::vector<std::vector<double>>(n + 1, std::vector<double>(max_span + 1, 0.0));

		for (int i = 1; i <= n; ++i)for (int j = i; j <= n && j - i <= max_span; ++j) {
			const WideFloating tmp = (GetZb(i, j) * GetWb(i, j));
			bppm[i][j - i] = exp(tmp.log_scale - total.log_scale);
		}
		std::vector<std::vector<int>>pseudo_bppm(n + 1, std::vector<int>(max_span + 1, 0));

		constexpr int size__ = 100000;
		for (int x = 0; x < size__; x++) {
			const std::string structure = GetRandomStructure();
			const auto p = VerifyAndParseStructure(structure, sequence, max_span, max_loop);
			for (int i = 1; i <= n; ++i)for (int j = i + 1; j <= n && j - i <= max_span; ++j)pseudo_bppm[i][j - i] += p[i][j];
		}
		for (int i = 1; i <= n; ++i)for (int j = i + 1; j <= n && j - i <= max_span; ++j) {

			const double diff = abs(double(pseudo_bppm[i][j - i]) / double(size__) - bppm[i][j - i]);
			if (diff >= exp(total.log_scale) * 0.01) {
				std::cout << std::endl;
				std::cout << sequence << std::endl;
				std::cout << i << std::endl;
				std::cout << j << std::endl;
				std::cout << max_loop << std::endl;
				std::cout << diff << std::endl;
				std::cout << pseudo_bppm[i][j - i] << std::endl;
				std::cout << total.log_scale << std::endl;
				std::cout << bppm[i][j - i] << std::endl;
				assert(0);
			}
		}
	}

	std::vector<std::string>ans;

	for (int i = 0; i < sample_amount; i++) {
		const std::string structure = GetRandomStructure();
		VerifyAndParseStructure(structure, sequence, max_span, max_loop);
		ans.push_back(structure);
	}

	return std::make_pair(ans, total);

}

}
