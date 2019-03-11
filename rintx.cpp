/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#include"rintx.h"

#include"parameter.h"
#include"misc.h"
#include"fourier_transform.h"
#include"complex_number.h"
#include"real_logsumexp.h"

namespace rintdwr {

//explicit instantiation

template std::vector<UsualComplexNumber<IntervalVar>> ComputeRintD1Dim(const RintX1DOptions& options);
template std::vector<UsualComplexNumber<Floating>> ComputeRintD1Dim(const RintX1DOptions& options);
template std::vector<WideComplexNumber<IntervalVar>> ComputeRintD1Dim(const RintX1DOptions& options);
template std::vector<WideComplexNumber<Floating>> ComputeRintD1Dim(const RintX1DOptions& options);

template std::vector<std::vector<UsualComplexNumber<IntervalVar>>>ComputeRintD2Dim(const RintX2DOptions& options);
template std::vector<std::vector<UsualComplexNumber<Floating>>>ComputeRintD2Dim(const RintX2DOptions& options);
template std::vector<std::vector<WideComplexNumber<IntervalVar>>>ComputeRintD2Dim(const RintX2DOptions& options);
template std::vector<std::vector<WideComplexNumber<Floating>>>ComputeRintD2Dim(const RintX2DOptions& options);

template std::pair<std::vector<UsualComplexNumber<IntervalVar>>, std::vector<std::vector<std::vector<UsualComplexNumber<IntervalVar>>>>> ComputeRintW1Dim(const RintX1DOptions& options);
template std::pair<std::vector<UsualComplexNumber<Floating>>, std::vector<std::vector<std::vector<UsualComplexNumber<Floating>>>>> ComputeRintW1Dim(const RintX1DOptions& options);
template std::pair<std::vector<WideComplexNumber<IntervalVar>>, std::vector<std::vector<std::vector<WideComplexNumber<IntervalVar>>>>> ComputeRintW1Dim(const RintX1DOptions& options);
template std::pair<std::vector<WideComplexNumber<Floating>>, std::vector<std::vector<std::vector<WideComplexNumber<Floating>>>>> ComputeRintW1Dim(const RintX1DOptions& options);

template std::pair<std::vector<UsualComplexNumber<IntervalVar>>, std::vector<std::vector<std::vector<UsualComplexNumber<IntervalVar>>>>> ComputeRintW1DimNonFourier(const RintX1DOptions& options);
template std::pair<std::vector<UsualComplexNumber<Floating>>, std::vector<std::vector<std::vector<UsualComplexNumber<Floating>>>>> ComputeRintW1DimNonFourier(const RintX1DOptions& options);
template std::pair<std::vector<WideComplexNumber<IntervalVar>>, std::vector<std::vector<std::vector<WideComplexNumber<IntervalVar>>>>> ComputeRintW1DimNonFourier(const RintX1DOptions& options);
template std::pair<std::vector<WideComplexNumber<Floating>>, std::vector<std::vector<std::vector<WideComplexNumber<Floating>>>>> ComputeRintW1DimNonFourier(const RintX1DOptions& options);

template std::pair<std::vector<std::vector<UsualComplexNumber<IntervalVar>>>, std::vector<std::vector<std::vector<std::vector<UsualComplexNumber<IntervalVar>>>>>> ComputeRintW2Dim(const RintX2DOptions& options);
template std::pair<std::vector<std::vector<UsualComplexNumber<Floating>>>, std::vector<std::vector<std::vector<std::vector<UsualComplexNumber<Floating>>>>>> ComputeRintW2Dim(const RintX2DOptions& options);
template std::pair<std::vector<std::vector<WideComplexNumber<IntervalVar>>>, std::vector<std::vector<std::vector<std::vector<WideComplexNumber<IntervalVar>>>>>> ComputeRintW2Dim(const RintX2DOptions& options);
template std::pair<std::vector<std::vector<WideComplexNumber<Floating>>>, std::vector<std::vector<std::vector<std::vector<WideComplexNumber<Floating>>>>>> ComputeRintW2Dim(const RintX2DOptions& options);

//function definition

template<typename Comp>std::vector<Comp>ComputeRintD1Dim(const RintX1DOptions& options) {

	typedef decltype(Comp::real) RealScalar;

	parasor_param::InitializeParameter(options.param_file_name, options.temperature);

	const int max_dim_ = options.allow_fft ? Ceiling2Power(options.max_dim1 + 1) : (options.max_dim1 + 1);

	const int n = int(options.sequence.size());
	const RealScalar pi = acos(RealScalar(-1.0));

	//[Mori et al., 2014].supp.(S10)
	const std::vector<std::vector<int>>C = ComputePredistanceMatrix(options.S1);

	//[Mori et al., 2014].supp.{(S9),(S14)...S1(21)}
	const auto g1 = [&](const int i, const int j) {return C[i][j]; };
	const auto g2 = [&](const int i, const int j, const int k) {return C[i][j] - C[i][k] - C[k + 1][j]; };
	const auto g3 = [&](const int i, const int j, const int k) {return C[i][j] - C[i][k]; };
	const auto g4 = [&](const int i, const int j) {return C[i][j] + 1 - 2 * options.S1[i][j]; };
	const auto g5 = [&](const int i, const int j, const int k, const int l) {return C[i][j] - C[k][l] + 1 - 2 * options.S1[i][j]; };
	const auto g6 = [&](const int i, const int j, const int k) {return C[i][j] - C[i + 1][k - 1] - C[k][j - 1] + 1 - 2 * options.S1[i][j]; };
	const auto g7 = [&](const int i, const int j, const int k) {return C[i][j] - C[k][j]; };
	const auto g8 = [&](const int i, const int j, const int k) {return C[i][j] - C[i][k - 1] - C[k][j]; };
	const auto g9 = g3;

	const auto at = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n && j - i <= options.max_span);
		return i * (options.max_span + 1) + (j - i);
	};

	//[Mori et al., 2014].(6)
	std::vector<Comp>zeta(max_dim_, Comp(0.0, 0.0));

	//[Mori et al., 2014].{(19)...(23)}
	//x is one of the (max_dim_)-th root of 1.
#pragma omp parallel for schedule(dynamic, 1)
	for (int x = 0; x < max_dim_; ++x) {

//#pragma omp critical
//		{
//			std::cout << x << "/" << max_dim_ << std::endl;
//		}

		std::vector<Comp> root_of_unity(max_dim_);
		for (int n1 = 0; n1 < max_dim_; ++n1) {
			const int d = (x * n1) % max_dim_;
			root_of_unity[n1] = Comp::GetPolar(
				RealScalar(1.0),
				RealScalar(2.0) * pi * RealScalar(double(d)) / RealScalar(double(max_dim_)));
		}

		const auto GetRoot = [&](const int i) {
			assert(0 <= i && i < max_dim_);
			return root_of_unity[i];
		};

		//data structure for memoization.
		//First is false at initialization, when it is calculated, it is set to true and a value is entered for second.
		typedef std::pair<bool, Comp>MemComp;

		//[Mori et al., 2014].{(19)...(23)}.  1-origin
		std::vector<MemComp>Z(n + 1, std::make_pair(false, Comp(0.0, 0.0)));
		std::vector<MemComp>Z1((n + 1) * (options.max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));
		std::vector<MemComp>Zb((n + 1) * (options.max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));
		std::vector<MemComp>Zm((n + 1) * (options.max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));
		std::vector<MemComp>Zm1((n + 1) * (options.max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));

		std::function<Comp(int)> GetZ;
		std::function<Comp(int, int)> GetZ1;
		std::function<Comp(int, int)> GetZb;
		std::function<Comp(int, int)> GetZm;
		std::function<Comp(int, int)> GetZm1;

		GetZ = [&](const int j) {
			assert(0 <= j && j <= n);
			if (j <= 1)return Comp(1.0, 0.0);
			if (Z[j].first)return Z[j].second;

			//[Mori et al., 2014].(19) is wrong. That is not possible to count cases where (1, *) forms base pair.
			Comp ans = GetZ1(1, j);
			ans += GetRoot(g1(1, j));

			for (int k = 1; k <= j - 1; ++k)ans += GetZ(k) * GetZ1(k + 1, j) * GetRoot(g2(1, j, k));

			return (Z[j] = std::make_pair(true, ans)).second;
		};
		GetZ1 = [&](const int i, const int j) {
			assert(1 <= i && i <= j && j <= n);
			if (i == j)return Comp(0.0, 0.0);

			if ((j - i) > options.max_span) {
				return GetZ1(i, i + options.max_span) * GetRoot(g3(i, j, i + options.max_span));
			}

			if (Z1[at(i, j)].first)return Z1[at(i, j)].second;

			Comp ans(0.0, 0.0);
			const int type = parasor_param::GetPairType(options.sequence[i - 1], options.sequence[j - 1]);
			if (type != 0) {
				ans += GetZb(i, j)
					* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, true, options.sequence)));
			}
			ans += GetZ1(i, j - 1) * GetRoot(g3(i, j, j - 1));
			return (Z1[at(i, j)] = std::make_pair(true, ans)).second;
		};
		GetZb = [&](const int i, const int j) {
			assert(1 <= i && i <= j && j <= n);
			if (i == j)return Comp(0.0, 0.0);

			const int type = parasor_param::GetPairType(options.sequence[i - 1], options.sequence[j - 1]);
			if (type == 0 || i + TURN >= j || (j - i) > options.max_span) {
				return Comp(0.0, 0.0);
			}

			if (Zb[at(i, j)].first)return Zb[at(i, j)].second;

			Comp ans(0.0, 0.0);

			//hairpin loop
			ans += Comp::LogRealToComp(RealScalar(parasor_param::ParHairpinEnergy(i - 1, j - 1, options.sequence))) * GetRoot(g4(i, j));

			//internal loop, stem, bulge
			for (int k = i + 1; k <= std::min(i + options.max_loop + 1, j - TURN - 2); ++k) {
				const int unpaired_base1 = k - i - 1;
				for (int l = std::max(k + TURN + 1, j - 1 - options.max_loop + unpaired_base1); l < j; ++l) {
					const int type = parasor_param::GetPairType(options.sequence[k - 1], options.sequence[l - 1]);
					if (type == 0)continue;
					ans += GetZb(k, l)
						* Comp::LogRealToComp(RealScalar(parasor_param::ParLoopEnergy(i - 1, j - 1, k - 1, l - 1, options.sequence)))
						* GetRoot(g5(i, j, k, l));
				}
			}

			//multi loop
			const int rtype = parasor_param::GetPairTypeReverse(options.sequence[i - 1], options.sequence[j - 1]);
			for (int k = i + TURN + 3; k + TURN + 1 <= j - 1; ++k) {
				ans += GetZm(i + 1, k - 1)
					* GetZm1(k, j - 1)
					* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(rtype, (j - 1) - 1, (i - 1) + 1, false, options.sequence)))
					* Comp::LogRealToComp(RealScalar(parasor_param::ParMultiloopClosing()))
					* GetRoot(g6(i, j, k));
			}

			return (Zb[at(i, j)] = std::make_pair(true, ans)).second;
		};
		GetZm = [&](const int i, const int j) {
			assert(1 <= i && (i <= j || j == i - 1) && j <= n);
			if (i == j)return Comp(0.0, 0.0);
			if (j == i - 1)return Comp(0.0, 0.0);
			if (Zm[at(i, j)].first)return Zm[at(i, j)].second;

			Comp ans(0.0, 0.0);
			for (int k = i; k + TURN + 1 <= j; ++k) {

				ans += (GetRoot(g7(i, j, k)) + GetZm(i, k - 1) * GetRoot(g8(i, j, k))) * GetZm1(k, j);
				//In [Mori et al., 2014], g7*exp(f4).
				//f4 is omitted because unpaired base is not used for energy calculation.

			}
			return (Zm[at(i, j)] = std::make_pair(true, ans)).second;
		};
		GetZm1 = [&](const int i, const int j) {
			assert(1 <= i && i <= j && j <= n);
			if (i == j)return Comp(0.0, 0.0);
			if (Zm1[at(i, j)].first)return Zm1[at(i, j)].second;

			Comp ans(0.0, 0.0);
			const int type = parasor_param::GetPairType(options.sequence[i - 1], options.sequence[j - 1]);
			if (type != 0) {
				ans += GetZb(i, j)
					* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, false, options.sequence)))
					* Comp::LogRealToComp(RealScalar(parasor_param::ParMultiloopInternal()));
			}
			ans += GetZm1(i, j - 1) * GetRoot(g9(i, j, j - 1));
			return (Zm1[at(i, j)] = std::make_pair(true, ans)).second;
		};

		zeta[x] = GetZ(n);
	}

	//[Mori et al., 2014].(6)
	const auto z = FourierTransform(zeta, options.allow_fft);

	if (options.allow_fft && max_dim_ != options.max_dim1) {
		std::vector<Comp>zz(options.max_dim1 + 1, Comp(0.0, 0.0));
		for (int d = 0; d <= options.max_dim1; ++d)zz[d] = z[d];
		return zz;
	}

	return z;
}
template<typename Comp>std::vector<std::vector<Comp>>ComputeRintD2Dim(const RintX2DOptions& options) {

	typedef decltype(Comp::real) RealScalar;

	parasor_param::InitializeParameter(options.param_file_name, options.temperature);

	const int max_dim_1 = options.allow_fft ? Ceiling2Power(options.max_dim1 + 1) : (options.max_dim1 + 1);
	const int max_dim_2 = options.allow_fft ? Ceiling2Power(options.max_dim2 + 1) : (options.max_dim2 + 1);

	const int n = int(options.sequence.size());
	const RealScalar pi = acos(RealScalar(-1.0));

	const std::vector<std::vector<int>> C1 = ComputePredistanceMatrix(options.S1);
	const std::vector<std::vector<int>> C2 = ComputePredistanceMatrix(options.S2);

	const auto g1_1 = [&](const int i, const int j) {return C1[i][j]; };
	const auto g2_1 = [&](const int i, const int j, const int k) {return C1[i][j] - C1[i][k] - C1[k + 1][j]; };
	const auto g3_1 = [&](const int i, const int j, const int k) {return C1[i][j] - C1[i][k]; };
	const auto g4_1 = [&](const int i, const int j) {return C1[i][j] + 1 - 2 * options.S1[i][j]; };
	const auto g5_1 = [&](const int i, const int j, const int k, const int l) {return C1[i][j] - C1[k][l] + 1 - 2 * options.S1[i][j]; };
	const auto g6_1 = [&](const int i, const int j, const int k) {return C1[i][j] - C1[i + 1][k - 1] - C1[k][j - 1] + 1 - 2 * options.S1[i][j]; };
	const auto g7_1 = [&](const int i, const int j, const int k) {return C1[i][j] - C1[k][j]; };
	const auto g8_1 = [&](const int i, const int j, const int k) {return C1[i][j] - C1[i][k - 1] - C1[k][j]; };
	const auto g9_1 = g3_1;

	const auto g1_2 = [&](const int i, const int j) {return C2[i][j]; };
	const auto g2_2 = [&](const int i, const int j, const int k) {return C2[i][j] - C2[i][k] - C2[k + 1][j]; };
	const auto g3_2 = [&](const int i, const int j, const int k) {return C2[i][j] - C2[i][k]; };
	const auto g4_2 = [&](const int i, const int j) {return C2[i][j] + 1 - 2 * options.S2[i][j]; };
	const auto g5_2 = [&](const int i, const int j, const int k, const int l) {return C2[i][j] - C2[k][l] + 1 - 2 * options.S2[i][j]; };
	const auto g6_2 = [&](const int i, const int j, const int k) {return C2[i][j] - C2[i + 1][k - 1] - C2[k][j - 1] + 1 - 2 * options.S2[i][j]; };
	const auto g7_2 = [&](const int i, const int j, const int k) {return C2[i][j] - C2[k][j]; };
	const auto g8_2 = [&](const int i, const int j, const int k) {return C2[i][j] - C2[i][k - 1] - C2[k][j]; };
	const auto g9_2 = g3_2;

	const auto at = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n && j - i <= options.max_span);
		return i * (options.max_span + 1) + (j - i);
	};

	std::vector<std::vector<Comp>>zeta(max_dim_1, std::vector<Comp>(max_dim_2, Comp(0.0, 0.0)));

#pragma omp parallel for schedule(dynamic, 1)
	for (int x = 0; x < max_dim_1 * max_dim_2; ++x) {
		const int x1 = x % max_dim_1;
		const int x2 = x / max_dim_1;

		std::vector<std::vector<Comp>>root_of_unity(max_dim_1, std::vector<Comp>(max_dim_2, Comp(0.0, 0.0)));
		for (int n1 = 0; n1 < max_dim_1; ++n1) {
			for (int n2 = 0; n2 < max_dim_2; ++n2) {
				const int d1 = (x1 * n1) % max_dim_1;
				const int d2 = (x2 * n2) % max_dim_2;
				root_of_unity[n1][n2] = Comp::GetPolar(
					RealScalar(1.0),
					RealScalar(2.0) * pi * (RealScalar(double(d1)) / RealScalar(double(max_dim_1)) + RealScalar(double(d2)) / RealScalar(double(max_dim_2))));
			}
		}

		const auto GetRoot = [&](const int i, const int j) {
			assert(0 <= i && i < max_dim_1 && 0 <= j && j < max_dim_2);
			return root_of_unity[i][j];
		};

		typedef std::pair<bool, Comp>MemComp;

		std::vector<MemComp>Z(n + 1, std::make_pair(false, Comp(0.0, 0.0)));
		std::vector<MemComp>Z1((n + 1) * (options.max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));
		std::vector<MemComp>Zb((n + 1) * (options.max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));
		std::vector<MemComp>Zm((n + 1) * (options.max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));
		std::vector<MemComp>Zm1((n + 1) * (options.max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));

		std::function<Comp(int)> GetZ;
		std::function<Comp(int, int)> GetZ1;
		std::function<Comp(int, int)> GetZb;
		std::function<Comp(int, int)> GetZm;
		std::function<Comp(int, int)> GetZm1;

		GetZ = [&](const int j) {
			assert(0 <= j && j <= n);
			if (j <= 1)return Comp(1.0, 0.0);
			if (Z[j].first)return Z[j].second;

			Comp ans = GetZ1(1, j);
			ans += GetRoot(g1_1(1, j), g1_2(1, j));

			for (int k = 1; k <= j - 1; ++k) {
				ans += GetZ(k) * GetZ1(k + 1, j) * GetRoot(g2_1(1, j, k), g2_2(1, j, k));
			}

			return (Z[j] = std::make_pair(true, ans)).second;
		};
		GetZ1 = [&](const int i, const int j) {
			assert(1 <= i && i <= j && j <= n);
			if (i == j)return Comp(0.0, 0.0);

			if ((j - i) > options.max_span) {
				return GetZ1(i, i + options.max_span) * GetRoot(g3_1(i, j, i + options.max_span), g3_2(i, j, i + options.max_span));
			}

			if (Z1[at(i, j)].first)return Z1[at(i, j)].second;

			Comp ans(0.0, 0.0);
			const int type = parasor_param::GetPairType(options.sequence[i - 1], options.sequence[j - 1]);
			if (type != 0) {
				ans += GetZb(i, j)
					* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, true, options.sequence)));
			}
			ans += GetZ1(i, j - 1) * GetRoot(g3_1(i, j, j - 1), g3_2(i, j, j - 1));

			return (Z1[at(i, j)] = std::make_pair(true, ans)).second;
		};
		GetZb = [&](const int i, const int j) {
			assert(1 <= i && i <= j && j <= n);
			if (i == j)return Comp(0.0, 0.0);
			if (Zb[at(i, j)].first)return Zb[at(i, j)].second;

			const int type = parasor_param::GetPairType(options.sequence[i - 1], options.sequence[j - 1]);
			if (type == 0 || i + TURN >= j || (j - i) > options.max_span) {
				return  Comp(0.0, 0.0);
			}

			Comp ans(0.0, 0.0);

			//hairpin loop
			ans += Comp::LogRealToComp(RealScalar(parasor_param::ParHairpinEnergy(i - 1, j - 1, options.sequence)))
				* GetRoot(g4_1(i, j), g4_2(i, j));

			//internal loop
			for (int k = i + 1; k <= std::min(i + options.max_loop + 1, j - TURN - 2); ++k) {
				const int unpaired_base1 = k - i - 1;
				for (int l = std::max(k + TURN + 1, j - 1 - options.max_loop + unpaired_base1); l < j; ++l) {
					const int type = parasor_param::GetPairType(options.sequence[k - 1], options.sequence[l - 1]);
					if (type == 0)continue;
					ans += GetZb(k, l)
						* Comp::LogRealToComp(RealScalar(parasor_param::ParLoopEnergy(i - 1, j - 1, k - 1, l - 1, options.sequence)))
						* GetRoot(g5_1(i, j, k, l), g5_2(i, j, k, l));
				}
			}

			//multi loop
			const int rtype = parasor_param::GetPairTypeReverse(options.sequence[i - 1], options.sequence[j - 1]);
			for (int k = i + TURN + 3; k + TURN + 1 <= j - 1; ++k) {
				ans += GetZm(i + 1, k - 1)
					* GetZm1(k, j - 1)
					* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(rtype, (j - 1) - 1, (i - 1) + 1, false, options.sequence)))
					* Comp::LogRealToComp(RealScalar(parasor_param::ParMultiloopClosing()))
					* GetRoot(g6_1(i, j, k), g6_2(i, j, k));
			}

			return (Zb[at(i, j)] = std::make_pair(true, ans)).second;
		};
		GetZm = [&](const int i, const int j) {
			assert(1 <= i && (i <= j || j == i - 1) && j <= n);
			if (i == j)return Comp(0.0, 0.0);
			if (j == i - 1)return Comp(0.0, 0.0);
			if (Zm[at(i, j)].first)return Zm[at(i, j)].second;

			Comp ans(0.0, 0.0);
			for (int k = i; k + TURN + 1 <= j; ++k) {
				ans += (/*Comp(1.0, 0.0) * */GetRoot(g7_1(i, j, k), g7_2(i, j, k)) + GetZm(i, k - 1) * GetRoot(g8_1(i, j, k), g8_2(i, j, k))) * GetZm1(k, j);
			}
			return (Zm[at(i, j)] = std::make_pair(true, ans)).second;
		};
		GetZm1 = [&](const int i, const int j) {
			assert(1 <= i && i <= j && j <= n);
			if (i == j)return Comp(0.0, 0.0);
			if (Zm1[at(i, j)].first)return Zm1[at(i, j)].second;

			Comp ans(0.0, 0.0);
			const int type = parasor_param::GetPairType(options.sequence[i - 1], options.sequence[j - 1]);
			if (type != 0) {
				ans += GetZb(i, j)
					* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, false, options.sequence)))
					* Comp::LogRealToComp(RealScalar(parasor_param::ParMultiloopInternal()));
			}
			ans += GetZm1(i, j - 1) * GetRoot(g9_1(i, j, j - 1), g9_2(i, j, j - 1));

			return (Zm1[at(i, j)] = std::make_pair(true, ans)).second;
		};

		zeta[x1][x2] = GetZ(n);
	}

	//hack: The following part is a little spaguetti.

	std::vector<std::vector<Comp>>z_dash(max_dim_1, std::vector<Comp>(max_dim_2, Comp(0.0, 0.0)));
	for (int x1 = 0; x1 < max_dim_1; ++x1) {
		std::vector<Comp>z_temp_a(max_dim_2, Comp(0.0, 0.0));
		for (int x2 = 0; x2 < max_dim_2; ++x2)z_temp_a[x2] = zeta[x1][x2];
		const auto z_temp_b = FourierTransform(z_temp_a, options.allow_fft);
		for (int x2 = 0; x2 < max_dim_2; ++x2) {
			z_dash[x1][x2] = z_temp_b[x2];
		}
	}

	std::vector<std::vector<Comp>>z(max_dim_1, std::vector<Comp>(max_dim_2, Comp(0.0, 0.0)));
	for (int x2 = 0; x2 < max_dim_2; ++x2) {
		std::vector<Comp>z_temp_a(max_dim_1, Comp(0.0, 0.0));
		for (int x1 = 0; x1 < max_dim_1; ++x1)z_temp_a[x1] = z_dash[x1][x2];
		const auto z_temp_b = FourierTransform(z_temp_a, options.allow_fft);
		for (int x1 = 0; x1 < max_dim_1; ++x1) {
			z[x1][x2] = z_temp_b[x1];
		}
	}

	if (options.allow_fft && (max_dim_1 != options.max_dim1 || max_dim_2 != options.max_dim2)) {
		std::vector<std::vector<Comp>>zz(options.max_dim1 + 1, std::vector<Comp>(options.max_dim2 + 1, Comp(0.0, 0.0)));
		for (int d1 = 0; d1 <= options.max_dim1; ++d1)for (int d2 = 0; d2 <= options.max_dim2; ++d2)zz[d1][d2] = z[d1][d2];
		return zz;
	}


	return z;
}
template<typename Comp>std::pair<std::vector<Comp>, std::vector<std::vector<std::vector<Comp>>>> ComputeRintW1Dim(const RintX1DOptions& options) {

	typedef decltype(Comp::real) RealScalar;

	parasor_param::InitializeParameter(options.param_file_name, options.temperature);

	const int max_dim_ = options.allow_fft ? Ceiling2Power(options.max_dim1 + 1) : (options.max_dim1 + 1);

	const int n = int(options.sequence.size());
	const RealScalar pi = acos(RealScalar(-1.0));

	const std::vector<std::vector<int>>C = ComputePredistanceMatrix(options.S1);

	const auto gz1 = [&](const int i, const int j) {return C[i][j]; };
	const auto gz2 = [&](const int i, const int j, const int k) {return C[i][j] - C[i][k] - C[k + 1][j]; };
	const auto gz3 = [&](const int i, const int j, const int k) {return C[i][j] - C[i][k]; };
	const auto gz4 = [&](const int i, const int j) {return C[i][j] + 1 - 2 * options.S1[i][j]; };
	const auto gz5 = [&](const int i, const int j, const int k, const int l) {return C[i][j] - C[k][l] + 1 - 2 * options.S1[i][j]; };
	const auto gz6 = [&](const int i, const int j, const int k) {return C[i][j] - C[i + 1][k - 1] - C[k][j - 1] + 1 - 2 * options.S1[i][j]; };
	const auto gz7 = [&](const int i, const int j, const int k) {return C[i][j] - C[k][j]; };
	const auto gz8 = [&](const int i, const int j, const int k) {return C[i][j] - C[i][k - 1] - C[k][j]; };
	const auto gz9 = gz3;

	//[Hagio et al., 2018].(4)
	//Since there is a possibility that i == n + 1, an "if" statement is added.
	const auto gz0 = [&](const int i, const int j) {
		if (i > j)return 0;
		return C[i][j];
	};

	//[Hagio et al., 2018].(14)
	const auto gw1 = [&](const int i, const int j) {return gz0(1, n) - gz0(i, j) - gz0(1, i - 1) - gz0(j + 1, n); };
	const auto gw2 = [&](const int i, const int j, const int h, const int l) {return gz0(h, l) - gz0(i, j) + 1 - 2 * options.S1[h][l]; };
	const auto gw3 = [&](const int i, const int j, const int h, const int l) {return gw2(i, j, h, l) - gz0(h + 1, i - 1); };
	const auto gw4 = [&](const int i, const int j, const int h, const int l) {return gw2(i, j, h, l) - gz0(j + 1, l - 1); };
	const auto gw5 = [&](const int i, const int j, const int h, const int l) {return gw2(i, j, h, l) - gz0(h + 1, i - 1) - gz0(j + 1, l - 1); };

	const auto at = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n && j - i <= options.max_span);
		return i * (options.max_span + 1) + (j - i);
	};

	std::vector<Comp>zeta(max_dim_, Comp(0.0, 0.0));
	std::vector<std::vector<Comp>>prob((n + 1) * (options.max_span + 1), std::vector<Comp>(max_dim_, Comp(0.0, 0.0)));

#pragma omp parallel for schedule(dynamic, 1)
	for (int x = 0; x < max_dim_; ++x) {
		
//#pragma omp critical
//		{
//			std::cout << x << "/" << max_dim_ << std::endl;
//		}

		std::vector<Comp> root_of_unity(max_dim_);
		for (int n1 = 0; n1 < max_dim_; ++n1) {
			const int d = (x * n1) % max_dim_;
			root_of_unity[n1] = Comp::GetPolar(
				RealScalar(1.0),
				RealScalar(2.0) * pi * RealScalar(double(d)) / RealScalar(double(max_dim_)));
		}

		const auto GetRoot = [&](const int i) {
			assert(0 <= i && i < max_dim_);
			return root_of_unity[i];
		};

		typedef std::pair<bool, Comp>MemComp;

		std::vector<MemComp>Z((n + 1) * 2, std::make_pair(false, Comp(0.0, 0.0)));
		std::vector<MemComp>Z1((n + 1) * (options.max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));
		std::vector<MemComp>Zb((n + 1) * (options.max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));
		std::vector<MemComp>Zm((n + 1) * (options.max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));
		std::vector<MemComp>Zm1((n + 1) * (options.max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));

		std::function<Comp(int, int)> GetZ;
		std::function<Comp(int, int)> GetZ1;
		std::function<Comp(int, int)> GetZb;
		std::function<Comp(int, int)> GetZm;
		std::function<Comp(int, int)> GetZm1;

		GetZ = [&](const int i, const int j) {
			assert(1 <= i && (i <= j || j == i - 1) && j <= n);
			if (i == j)return Comp(1.0, 0.0);
			if (j == i - 1)return Comp(1.0, 0.0);
			assert(i == 1 || j == n);
			const int index = (i == 1) ? j : (i + n + 1);
			assert(1 <= index && index < (n + 1) * 2);
			if (Z[index].first)return Z[index].second;
			Comp ans(0.0, 0.0);

			if (i == 1) {

				//no base pair
				ans += GetRoot(gz1(i, j));

				//there is only one outermost base pair (i, *)
				ans += GetZ1(i, j);

				//there are one or more outermost base pairs. The rightmost outermost base pair is (k + 1, *)
				for (int k = i; k <= j - 1; ++k) {
					ans += GetZ(i, k) * GetZ1(k + 1, j) * GetRoot(/*gz2(i, j, k)*/gz0(i, j) - gz0(i, k) - gz0(k + 1, j));
				}
			}
			else {

				//i is not paired
				ans += GetZ(i + 1, j) * GetRoot(gz0(i, j) - gz0(i + 1, j));

				//(i,k) form a base pair
				for (int k = i + TURN + 1; k <= j && k - i <= options.max_span; ++k) {
					const int type = parasor_param::GetPairType(options.sequence[i - 1], options.sequence[k - 1]);
					if (type == 0)continue;
					ans += GetZb(i, k)
						* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (k - 1) + 1, true, options.sequence)))
						* GetZ(k + 1, j)
						* GetRoot(/*gz2(i, j, k)*/gz0(i, j) - gz0(i, k) - gz0(k + 1, j));
				}
			}
			return (Z[index] = std::make_pair(true, ans)).second;
		};
		GetZ1 = [&](const int i, const int j) {
			assert(1 <= i && i <= j && j <= n);
			if (i == j)return Comp(0.0, 0.0);

			if ((j - i) > options.max_span) {
				return GetZ1(i, i + options.max_span) * GetRoot(gz3(i, j, i + options.max_span));
			}

			if (Z1[at(i, j)].first)return Z1[at(i, j)].second;

			Comp ans(0.0, 0.0);
			const int type = parasor_param::GetPairType(options.sequence[i - 1], options.sequence[j - 1]);
			if (type != 0) {
				ans += GetZb(i, j)
					* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, true, options.sequence)));
			}
			ans += GetZ1(i, j - 1) * GetRoot(gz3(i, j, j - 1));

			return (Z1[at(i, j)] = std::make_pair(true, ans)).second;
		};
		GetZb = [&](const int i, const int j) {
			assert(1 <= i && i <= j && j <= n);
			if (i == j)return Comp(0.0, 0.0);

			const int type = parasor_param::GetPairType(options.sequence[i - 1], options.sequence[j - 1]);
			if (type == 0 || i + TURN >= j || (j - i) > options.max_span) {
				return Comp(0.0, 0.0);
			}

			if (Zb[at(i, j)].first)return Zb[at(i, j)].second;

			Comp ans(0.0, 0.0);

			//hairpin loop
			ans += Comp::LogRealToComp(RealScalar(parasor_param::ParHairpinEnergy(i - 1, j - 1, options.sequence))) * GetRoot(gz4(i, j));

			//internal loop
			for (int k = i + 1; k <= std::min(i + options.max_loop + 1, j - TURN - 2); ++k) {
				const int unpaired_base1 = k - i - 1;
				for (int l = std::max(k + TURN + 1, j - 1 - options.max_loop + unpaired_base1); l < j; ++l) {
					const int type = parasor_param::GetPairType(options.sequence[k - 1], options.sequence[l - 1]);
					if (type == 0)continue;
					ans += GetZb(k, l)
						* Comp::LogRealToComp(RealScalar(parasor_param::ParLoopEnergy(i - 1, j - 1, k - 1, l - 1, options.sequence)))
						* GetRoot(gz5(i, j, k, l));
				}
			}

			//multi loop
			const int rtype = parasor_param::GetPairTypeReverse(options.sequence[i - 1], options.sequence[j - 1]);
			for (int k = i + TURN + 3; k + TURN + 1 <= j - 1; ++k) {
				ans += GetZm(i + 1, k - 1)
					* GetZm1(k, j - 1)
					* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(rtype, (j - 1) - 1, (i - 1) + 1, false, options.sequence)))
					* Comp::LogRealToComp(RealScalar(parasor_param::ParMultiloopClosing()))
					* GetRoot(gz6(i, j, k));
			}

			return (Zb[at(i, j)] = std::make_pair(true, ans)).second;
		};
		GetZm = [&](const int i, const int j) {
			assert(1 <= i && (i <= j || j == i - 1) && j <= n);
			if (i == j)return Comp(0.0, 0.0);
			if (j == i - 1)return Comp(0.0, 0.0);
			if (Zm[at(i, j)].first)return Zm[at(i, j)].second;

			Comp ans(0.0, 0.0);
			for (int k = i; k + TURN + 1 <= j; ++k) {
				ans += (GetRoot(gz7(i, j, k)) + GetZm(i, k - 1) * GetRoot(gz8(i, j, k))) * GetZm1(k, j);
			}
			return (Zm[at(i, j)] = std::make_pair(true, ans)).second;
		};
		GetZm1 = [&](const int i, const int j) {
			assert(1 <= i && i <= j && j <= n);
			if (i == j)return Comp(0.0, 0.0);
			if (Zm1[at(i, j)].first)return Zm1[at(i, j)].second;

			Comp ans(0.0, 0.0);
			const int type = parasor_param::GetPairType(options.sequence[i - 1], options.sequence[j - 1]);
			if (type != 0) {
				ans += GetZb(i, j)
					* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, false, options.sequence)))
					* Comp::LogRealToComp(RealScalar(parasor_param::ParMultiloopInternal()));
			}
			ans += GetZm1(i, j - 1) * GetRoot(gz9(i, j, j - 1));

			return (Zm1[at(i, j)] = std::make_pair(true, ans)).second;
		};

		//[Hagio et al., 2018], Algorithm2
		std::vector<MemComp>Wb((n + 1) * (options.max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));
		std::function<Comp(int, int)> GetWb;
		GetWb = [&](const int i, const int j) {
			assert(1 <= i && i <= j && j <= n);

			const int type = parasor_param::GetPairType(options.sequence[i - 1], options.sequence[j - 1]);
			if (type == 0 || i + TURN >= j || (j - i) > options.max_span) {
				return Comp(0.0, 0.0);
			}

			if (Wb[at(i, j)].first)return Wb[at(i, j)].second;

			Comp ans(0.0, 0.0);

			//(i,j) is an outermost base pair
			ans += GetZ(1, i - 1)
				* GetZ(j + 1, n)
				* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, true, options.sequence)))
				* GetRoot(gw1(i, j));

			//base pair (h,l) exists outside (i,j)
			for (int h = std::max(1, std::max(i - options.max_span + TURN + 2, i - options.max_loop - 1)); h < i; ++h) {
				for (int l = j + 1; l <= n && (l - h) <= options.max_span && (i - h - 1) + (l - j - 1) <= options.max_loop; ++l) {
					//internal loop, bulge, stem
					const int type = parasor_param::GetPairType(options.sequence[h - 1], options.sequence[l - 1]);
					if (type == 0)continue;
					ans += GetWb(h, l)
						* Comp::LogRealToComp(RealScalar(parasor_param::ParLoopEnergy(h - 1, l - 1, i - 1, j - 1, options.sequence)))
						* GetRoot(gw2(i, j, h, l));
				}
			}
			for (int h = std::max(1, i - options.max_span + TURN + 2); h < i; ++h) {
				for (int l = j + 1; l <= n && (l - h) <= options.max_span; ++l) {
					//(h,l) is multiloop closing base pair. (i,j) is one of the branches
					const int rtype = parasor_param::GetPairTypeReverse(options.sequence[h - 1], options.sequence[l - 1]);
					if (rtype == 0)continue;
					ans += Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, false, options.sequence)))
						* Comp::LogRealToComp(RealScalar(parasor_param::ParMultiloopInternal()))
						* GetWb(h, l)
						* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(rtype, (l - 1) - 1, (h - 1) + 1, false, options.sequence)))
						* Comp::LogRealToComp(RealScalar(parasor_param::ParMultiloopClosing()))
						* (
							GetZm(h + 1, i - 1) * GetRoot(gw3(i, j, h, l))
							+ GetZm(j + 1, l - 1) * GetRoot(gw4(i, j, h, l))
							+ GetZm(h + 1, i - 1) * GetZm(j + 1, l - 1) * GetRoot(gw5(i, j, h, l))
							);
				}
			}
			return (Wb[at(i, j)] = std::make_pair(true, ans)).second;
		};

		zeta[x] = GetZ(1, n);

		for (int i = 1; i <= n; ++i)for (int j = i; j <= n && (j - i) <= options.max_span; ++j) {

			//[Hagio et al., 2018].(9)
			prob[at(i, j)][x] = GetZb(i, j) * GetWb(i, j);
		}
	}

	//std::vector<std::vector<Comp>>prob((n + 1) * (options.max_span + 1), std::vector<Comp>(max_dim_, Comp(0.0, 0.0)));
	std::vector<std::vector<std::vector<Comp>>>P(
		options.max_dim1 + 1, std::vector<std::vector<Comp>>(
			n + 1, std::vector<Comp>(
				options.max_span + 1, Comp(0.0, 0.0))));

	std::vector<std::pair<int, int>>index;
	for (int i = 1; i <= n; ++i) for (int j = i + 1; j <= n && (j - i) <= options.max_span; ++j)index.push_back(std::make_pair(i, j));

	int count = 0;
#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < int(index.size()); ++i) {

//#pragma omp critical
//		{
//			if (i % 100 == 0) std::cout << i << "/" << int(index.size()) << std::endl;
//		}

		const int a = index[i].first;
		const int b = index[i].second;
		const std::vector<Comp>p_ = FourierTransform(prob[at(a, b)], options.allow_fft);
		for (int d = 0; d <= options.max_dim1; ++d)P[d][a][b - a] = p_[d];
	}

	const std::vector<Comp>z = FourierTransform(zeta, options.allow_fft);

	if (options.allow_fft && max_dim_ != options.max_dim1) {
		std::vector<Comp>zz(options.max_dim1 + 1, Comp(0.0, 0.0));
		for (int d = 0; d <= options.max_dim1; ++d)zz[d] = z[d];
		return std::make_pair(zz, P);
	}

	return std::make_pair(z, P);
}
template<typename Comp>std::pair<std::vector<Comp>, std::vector<std::vector<std::vector<Comp>>>> ComputeRintW1DimNonFourier(const RintX1DOptions& options) {

	//Fourier transform is not used. Time complexity is O(NW^3Hmax^3)

	//フーリエ変換を使わない。計算量はNW^3Hmax^3

	typedef decltype(Comp::real) RealScalar;

	parasor_param::InitializeParameter(options.param_file_name, options.temperature);

	const int n = int(options.sequence.size());

	const std::vector<std::vector<int>>C = ComputePredistanceMatrix(options.S1);

	const auto gz1 = [&](const int i, const int j) {return C[i][j]; };
	const auto gz2 = [&](const int i, const int j, const int k) {return C[i][j] - C[i][k] - C[k + 1][j]; };
	const auto gz3 = [&](const int i, const int j, const int k) {return C[i][j] - C[i][k]; };
	const auto gz4 = [&](const int i, const int j) {return C[i][j] + 1 - 2 * options.S1[i][j]; };
	const auto gz5 = [&](const int i, const int j, const int k, const int l) {return C[i][j] - C[k][l] + 1 - 2 * options.S1[i][j]; };
	const auto gz6 = [&](const int i, const int j, const int k) {return C[i][j] - C[i + 1][k - 1] - C[k][j - 1] + 1 - 2 * options.S1[i][j]; };
	const auto gz7 = [&](const int i, const int j, const int k) {return C[i][j] - C[k][j]; };
	const auto gz8 = [&](const int i, const int j, const int k) {return C[i][j] - C[i][k - 1] - C[k][j]; };
	const auto gz9 = gz3;

	const auto gz0 = [&](const int i, const int j) {
		if (i > j)return 0;
		return C[i][j];
	};

	const auto gw1 = [&](const int i, const int j) {return gz0(1, n) - gz0(i, j) - gz0(1, i - 1) - gz0(j + 1, n); };
	const auto gw2 = [&](const int i, const int j, const int h, const int l) {return gz0(h, l) - gz0(i, j) + 1 - 2 * options.S1[h][l]; };
	const auto gw3 = [&](const int i, const int j, const int h, const int l) {return gw2(i, j, h, l) - gz0(h + 1, i - 1); };
	const auto gw4 = [&](const int i, const int j, const int h, const int l) {return gw2(i, j, h, l) - gz0(j + 1, l - 1); };
	const auto gw5 = [&](const int i, const int j, const int h, const int l) {return gw2(i, j, h, l) - gz0(h + 1, i - 1) - gz0(j + 1, l - 1); };

	const auto at = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n && j - i <= options.max_span);
		return i * (options.max_span + 1) + (j - i);
	};

	typedef std::pair<bool, Comp>MemReal;

	std::vector<std::vector<MemReal>>Z(options.max_dim1 + 1, std::vector<MemReal>((n + 1) * 2, std::make_pair(false, Comp(0.0, 0.0))));
	std::vector<std::vector<MemReal>>Z1(options.max_dim1 + 1, std::vector<MemReal>((n + 1) * (options.max_span + 1), std::make_pair(false, Comp(0.0, 0.0))));
	std::vector<std::vector<MemReal>>Zb(options.max_dim1 + 1, std::vector<MemReal>((n + 1) * (options.max_span + 1), std::make_pair(false, Comp(0.0, 0.0))));
	std::vector<std::vector<MemReal>>Zm(options.max_dim1 + 1, std::vector<MemReal>((n + 1) * (options.max_span + 1), std::make_pair(false, Comp(0.0, 0.0))));
	std::vector<std::vector<MemReal>>Zm1(options.max_dim1 + 1, std::vector<MemReal>((n + 1) * (options.max_span + 1), std::make_pair(false, Comp(0.0, 0.0))));
	std::vector<std::vector<MemReal>>Wb(options.max_dim1 + 1, std::vector<MemReal>((n + 1) * (options.max_span + 1), std::make_pair(false, Comp(0.0, 0.0))));

	std::function<Comp(int, int, int)> GetZ;
	std::function<Comp(int, int, int)> GetZ1;
	std::function<Comp(int, int, int)> GetZb;
	std::function<Comp(int, int, int)> GetZm;
	std::function<Comp(int, int, int)> GetZm1;
	std::function<Comp(int, int, int)> GetWb;

	GetZ = [&](const int i, const int j, const int x) {
		assert(1 <= i && (i <= j || j == i - 1) && j <= n);
		assert(0 <= x && x <= options.max_dim1);
		if (i == j)return Comp(x == 0 ? 1.0 : 0.0, 0.0);
		if (j == i - 1)return Comp(x == 0 ? 1.0 : 0.0, 0.0);
		assert(i == 1 || j == n);
		const int index = (i == 1) ? j : (i + n + 1);
		assert(1 <= index && index < (n + 1) * 2);
		if (Z[x][index].first)return Z[x][index].second;

		Comp ans(0.0, 0.0);

		if (i == 1) {

			if (x == gz1(1, j))ans += Comp(1.0, 0.0);

			for (int k = i - 1; k <= j - 1; ++k) {
				const int y = x - (gz0(i, j) - gz0(i, k) - gz0(k + 1, j));
				for (int z = 0; z <= y; ++z) {
					ans += GetZ(1, k, z) * GetZ1(k + 1, j, y - z);
				}
			}
		}
		else {

			const int y1 = x - (gz0(i, j) - gz0(i + 1, j));
			if (y1 >= 0)ans += GetZ(i + 1, j, y1);

			for (int k = i + TURN + 1; k <= j && k - i <= options.max_span; ++k) {
				const int type = parasor_param::GetPairType(options.sequence[i - 1], options.sequence[k - 1]);
				if (type == 0)continue;
				const int y2 = x - (gz0(i, j) - gz0(i, k) - gz0(k + 1, j));
				for (int z = 0; z <= y2; ++z) {
					ans += GetZb(i, k, z)
						* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (k - 1) + 1, true, options.sequence)))
						* GetZ(k + 1, j, y2 - z);
				}
			}
		}
		return (Z[x][index] = std::make_pair(true, ans)).second;
	};
	GetZ1 = [&](const int i, const int j, const int x) {
		assert(1 <= i && i <= j && j <= n);
		assert(0 <= x && x <= options.max_dim1);
		if (i == j)return Comp(0.0, 0.0);

		if ((j - i) > options.max_span) {
			const int y = x - gz3(i, j, i + options.max_span);
			return y >= 0 ? GetZ1(i, i + options.max_span, y) : Comp(0.0, 0.0);
		}

		if (Z1[x][at(i, j)].first)return Z1[x][at(i, j)].second;

		Comp ans(0.0, 0.0);

		const int type = parasor_param::GetPairType(options.sequence[i - 1], options.sequence[j - 1]);
		if (type != 0) {
			ans += GetZb(i, j, x)
				* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, true, options.sequence)));
		}
		const int y = x - gz3(i, j, j - 1);
		if (y >= 0)ans += GetZ1(i, j - 1, y);

		return (Z1[x][at(i, j)] = std::make_pair(true, ans)).second;
	};
	GetZb = [&](const int i, const int j, const int x) {
		assert(1 <= i && i <= j && j <= n);
		assert(0 <= x && x <= options.max_dim1);
		if (i == j)return Comp(0.0, 0.0);

		const int type = parasor_param::GetPairType(options.sequence[i - 1], options.sequence[j - 1]);
		if (type == 0 || i + TURN >= j || (j - i) > options.max_span) {
			return Comp(0.0, 0.0);
		}

		if (Zb[x][at(i, j)].first)return Zb[x][at(i, j)].second;

		Comp ans(0.0, 0.0);

		//hairpin loop
		if (x == gz4(i, j))ans += Comp::LogRealToComp(RealScalar(parasor_param::ParHairpinEnergy(i - 1, j - 1, options.sequence)));

		//internal loop
		for (int k = i + 1; k <= std::min(i + options.max_loop + 1, j - TURN - 2); ++k) {
			const int unpaired_base1 = k - i - 1;
			for (int l = std::max(k + TURN + 1, j - 1 - options.max_loop + unpaired_base1); l < j; ++l) {
				const int type = parasor_param::GetPairType(options.sequence[k - 1], options.sequence[l - 1]);
				if (type == 0)continue;

				const int y = x - gz5(i, j, k, l);
				if (y >= 0) {
					ans += GetZb(k, l, y)
						* Comp::LogRealToComp(RealScalar(parasor_param::ParLoopEnergy(i - 1, j - 1, k - 1, l - 1, options.sequence)));

				}
			}
		}

		//multi loop
		const int rtype = parasor_param::GetPairTypeReverse(options.sequence[i - 1], options.sequence[j - 1]);
		for (int k = i + TURN + 3; k + TURN + 1 <= j - 1; ++k) {
			const int y = x - gz6(i, j, k);
			for (int z = 0; z <= y; ++z) {
				ans += GetZm(i + 1, k - 1, z)
					* GetZm1(k, j - 1, y - z)
					* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(rtype, (j - 1) - 1, (i - 1) + 1, false, options.sequence)))
					* Comp::LogRealToComp(RealScalar(parasor_param::ParMultiloopClosing()));
			}
		}

		return (Zb[x][at(i, j)] = std::make_pair(true, ans)).second;
	};
	GetZm = [&](const int i, const int j, const int x) {
		assert(1 <= i && (i <= j || j == i - 1) && j <= n);
		assert(0 <= x && x <= options.max_dim1);
		if (i == j)return Comp(0.0, 0.0);
		if (j == i - 1)return Comp(0.0, 0.0);
		if (Zm[x][at(i, j)].first)return Zm[x][at(i, j)].second;

		Comp ans(0.0, 0.0);
		for (int k = i; k + TURN + 1 <= j; ++k) {
			const int y1 = x - gz7(i, j, k);
			if (y1 >= 0)ans += GetZm1(k, j, y1);
			const int y2 = x - gz8(i, j, k);
			for (int z = 0; z <= y2; ++z)ans += GetZm(i, k - 1, z) * GetZm1(k, j, y2 - z);
		}

		return (Zm[x][at(i, j)] = std::make_pair(true, ans)).second;
	};
	GetZm1 = [&](const int i, const int j, const int x) {
		assert(1 <= i && i <= j && j <= n);
		assert(0 <= x && x <= options.max_dim1);
		if (i == j)return Comp(0.0, 0.0);
		if (Zm1[x][at(i, j)].first)return Zm1[x][at(i, j)].second;

		Comp ans(0.0, 0.0);
		for (int k = i + TURN + 1; k <= j; ++k) {
			const int type = parasor_param::GetPairType(options.sequence[i - 1], options.sequence[k - 1]);
			if (type == 0)continue;
			const int y = x - gz9(i, j, k);
			if (y >= 0) {
				ans += GetZb(i, k, y)
					* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (k - 1) + 1, false, options.sequence)))
					* Comp::LogRealToComp(RealScalar(parasor_param::ParMultiloopInternal()));
			}
		}
		return (Zm1[x][at(i, j)] = std::make_pair(true, ans)).second;
	};
	GetWb = [&](const int i, const int j, const int x) {
		assert(1 <= i && i <= j && j <= n);
		assert(0 <= x && x <= options.max_dim1);

		const int type = parasor_param::GetPairType(options.sequence[i - 1], options.sequence[j - 1]);
		if (type == 0 || i + TURN >= j || (j - i) > options.max_span) {
			return Comp(0.0, 0.0);
		}

		if (Wb[x][at(i, j)].first)return Wb[x][at(i, j)].second;

		Comp ans(0.0, 0.0);

		const int y1 = x - gw1(i, j);
		for (int z = 0; z <= y1; ++z) {
			ans += GetZ(1, i - 1, z)
				* GetZ(j + 1, n, y1 - z)
				* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, true, options.sequence)));
		}

		for (int h = std::max(1, std::max(i - options.max_span + TURN + 2, i - options.max_loop - 1)); h < i; ++h) {
			for (int l = j + 1; l <= n && (l - h) <= options.max_span && (i - h - 1) + (l - j - 1) <= options.max_loop; ++l) {
				//internal loop, bulge, stem
				const int type = parasor_param::GetPairType(options.sequence[h - 1], options.sequence[l - 1]);
				if (type == 0)continue;
				const int y2 = x - gw2(i, j, h, l);
				if (y2 >= 0) {
					ans += GetWb(h, l, y2)
						* Comp::LogRealToComp(RealScalar(parasor_param::ParLoopEnergy(h - 1, l - 1, i - 1, j - 1, options.sequence)));
				}
			}
		}
		for (int h = std::max(1, i - options.max_span + TURN + 2); h < i; ++h) {
			for (int l = j + 1; l <= n && (l - h) <= options.max_span; ++l) {

				const int rtype = parasor_param::GetPairTypeReverse(options.sequence[h - 1], options.sequence[l - 1]);
				if (rtype == 0)continue;

				const int y3 = x - gw3(i, j, h, l);
				for (int z = 0; z <= y3; ++z) {
					ans += Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, false, options.sequence)))
						* Comp::LogRealToComp(RealScalar(parasor_param::ParMultiloopInternal()))
						* GetWb(h, l, z)
						* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(rtype, (l - 1) - 1, (h - 1) + 1, false, options.sequence)))
						* Comp::LogRealToComp(RealScalar(parasor_param::ParMultiloopClosing()))
						* GetZm(h + 1, i - 1, y3 - z);
				}

				const int y4 = x - gw4(i, j, h, l);
				for (int z = 0; z <= y4; ++z) {
					ans += Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, false, options.sequence)))
						* Comp::LogRealToComp(RealScalar(parasor_param::ParMultiloopInternal()))
						* GetWb(h, l, z)
						* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(rtype, (l - 1) - 1, (h - 1) + 1, false, options.sequence)))
						* Comp::LogRealToComp(RealScalar(parasor_param::ParMultiloopClosing()))
						* GetZm(j + 1, l - 1, y4 - z);
				}

				const int y5 = x - gw5(i, j, h, l);
				for (int z1 = 0; z1 <= y5; ++z1)for (int z2 = 0; z1 + z2 <= y5; ++z2) {
					ans += Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, false, options.sequence)))
						* Comp::LogRealToComp(RealScalar(parasor_param::ParMultiloopInternal()))
						* GetWb(h, l, z1)
						* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(rtype, (l - 1) - 1, (h - 1) + 1, false, options.sequence)))
						* Comp::LogRealToComp(RealScalar(parasor_param::ParMultiloopClosing()))
						* GetZm(h + 1, i - 1, z2)
						* GetZm(j + 1, l - 1, y5 - z1 - z2);
				}
			}
		}
		return (Wb[x][at(i, j)] = std::make_pair(true, ans)).second;
	};

	std::vector<Comp>ans_Z(options.max_dim1 + 1, Comp(0.0, 0.0));
	std::vector<std::vector<std::vector<Comp>>>ans_P(
		options.max_dim1 + 1, std::vector<std::vector<Comp>>(
			n + 1, std::vector<Comp>(
				options.max_span + 1, Comp(0.0, 0.0))));

	for (int x = 0; x <= options.max_dim1; ++x) {

		//std::cout << "NonFourier: " << x << "/" << options.max_dim1 + 1 << std::endl;

		ans_Z[x] = GetZ(1, n, x);

		for (int i = 1; i <= n; ++i)for (int j = i; j <= n && (j - i) <= options.max_span; ++j) {
			for (int k = 0; k <= x; ++k)ans_P[x][i][j - i] += GetZb(i, j, k) * GetWb(i, j, x - k);
		}
	}

	return std::make_pair(ans_Z, ans_P);
}
template<typename Comp>std::pair<std::vector<std::vector<Comp>>, std::vector<std::vector<std::vector<std::vector<Comp>>>>> ComputeRintW2Dim(const RintX2DOptions& options) {

	typedef decltype(Comp::real) RealScalar;

	parasor_param::InitializeParameter(options.param_file_name, options.temperature);

	const int max_dim_1 = options.allow_fft ? Ceiling2Power(options.max_dim1 + 1) : (options.max_dim1 + 1);
	const int max_dim_2 = options.allow_fft ? Ceiling2Power(options.max_dim2 + 1) : (options.max_dim2 + 1);

	const int n = int(options.sequence.size());
	const RealScalar pi = acos(RealScalar(-1.0));

	const std::vector<std::vector<int>>C1 = ComputePredistanceMatrix(options.S1);
	const std::vector<std::vector<int>>C2 = ComputePredistanceMatrix(options.S2);

	const auto gz1_1 = [&](const int i, const int j) {return C1[i][j]; };
	const auto gz2_1 = [&](const int i, const int j, const int k) {return C1[i][j] - C1[i][k] - C1[k + 1][j]; };
	const auto gz3_1 = [&](const int i, const int j, const int k) {return C1[i][j] - C1[i][k]; };
	const auto gz4_1 = [&](const int i, const int j) {return C1[i][j] + 1 - 2 * options.S1[i][j]; };
	const auto gz5_1 = [&](const int i, const int j, const int k, const int l) {return C1[i][j] - C1[k][l] + 1 - 2 * options.S1[i][j]; };
	const auto gz6_1 = [&](const int i, const int j, const int k) {return C1[i][j] - C1[i + 1][k - 1] - C1[k][j - 1] + 1 - 2 * options.S1[i][j]; };
	const auto gz7_1 = [&](const int i, const int j, const int k) {return C1[i][j] - C1[k][j]; };
	const auto gz8_1 = [&](const int i, const int j, const int k) {return C1[i][j] - C1[i][k - 1] - C1[k][j]; };
	const auto gz9_1 = gz3_1;

	const auto gz1_2 = [&](const int i, const int j) {return C2[i][j]; };
	const auto gz2_2 = [&](const int i, const int j, const int k) {return C2[i][j] - C2[i][k] - C2[k + 1][j]; };
	const auto gz3_2 = [&](const int i, const int j, const int k) {return C2[i][j] - C2[i][k]; };
	const auto gz4_2 = [&](const int i, const int j) {return C2[i][j] + 1 - 2 * options.S2[i][j]; };
	const auto gz5_2 = [&](const int i, const int j, const int k, const int l) {return C2[i][j] - C2[k][l] + 1 - 2 * options.S2[i][j]; };
	const auto gz6_2 = [&](const int i, const int j, const int k) {return C2[i][j] - C2[i + 1][k - 1] - C2[k][j - 1] + 1 - 2 * options.S2[i][j]; };
	const auto gz7_2 = [&](const int i, const int j, const int k) {return C2[i][j] - C2[k][j]; };
	const auto gz8_2 = [&](const int i, const int j, const int k) {return C2[i][j] - C2[i][k - 1] - C2[k][j]; };
	const auto gz9_2 = gz3_2;

	//[Hagio et al., 2018].(4)
	//Since there is a possibility that i == n + 1, an "if" statement is added.
	const auto gz0_1 = [&](const int i, const int j) {
		if (i > j)return 0;
		return C1[i][j];
	};
	const auto gz0_2 = [&](const int i, const int j) {
		if (i > j)return 0;
		return C2[i][j];
	};

	//[Hagio et al., 2018].(14)
	const auto gw1_1 = [&](const int i, const int j) {return gz0_1(1, n) - gz0_1(i, j) - gz0_1(1, i - 1) - gz0_1(j + 1, n); };
	const auto gw2_1 = [&](const int i, const int j, const int h, const int l) {return gz0_1(h, l) - gz0_1(i, j) + 1 - 2 * options.S1[h][l]; };
	const auto gw3_1 = [&](const int i, const int j, const int h, const int l) {return gw2_1(i, j, h, l) - gz0_1(h + 1, i - 1); };
	const auto gw4_1 = [&](const int i, const int j, const int h, const int l) {return gw2_1(i, j, h, l) - gz0_1(j + 1, l - 1); };
	const auto gw5_1 = [&](const int i, const int j, const int h, const int l) {return gw2_1(i, j, h, l) - gz0_1(h + 1, i - 1) - gz0_1(j + 1, l - 1); };

	const auto gw1_2 = [&](const int i, const int j) {return gz0_2(1, n) - gz0_2(i, j) - gz0_2(1, i - 1) - gz0_2(j + 1, n); };
	const auto gw2_2 = [&](const int i, const int j, const int h, const int l) {return gz0_2(h, l) - gz0_2(i, j) + 1 - 2 * options.S2[h][l]; };
	const auto gw3_2 = [&](const int i, const int j, const int h, const int l) {return gw2_2(i, j, h, l) - gz0_2(h + 1, i - 1); };
	const auto gw4_2 = [&](const int i, const int j, const int h, const int l) {return gw2_2(i, j, h, l) - gz0_2(j + 1, l - 1); };
	const auto gw5_2 = [&](const int i, const int j, const int h, const int l) {return gw2_2(i, j, h, l) - gz0_2(h + 1, i - 1) - gz0_2(j + 1, l - 1); };

	const auto at = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n && j - i <= options.max_span);
		return i * (options.max_span + 1) + (j - i);
	};

	std::vector<std::vector<Comp>>zeta(max_dim_1, std::vector<Comp>(max_dim_2, Comp(0.0, 0.0)));
	std::vector<std::vector<std::vector<Comp>>>prob(
		(n + 1) * (options.max_span + 1), std::vector<std::vector<Comp>>(
			max_dim_1, std::vector<Comp>(
				max_dim_2, Comp(0.0, 0.0))));

#pragma omp parallel for schedule(dynamic, 1)
	for (int x = 0; x < max_dim_1 * max_dim_2; ++x) {

		const int x1 = x % max_dim_1;
		const int x2 = x / max_dim_1;

		std::vector<std::vector<Comp>>root_of_unity(max_dim_1, std::vector<Comp>(max_dim_2, Comp(0.0, 0.0)));
		for (int n1 = 0; n1 < max_dim_1; ++n1) {
			for (int n2 = 0; n2 < max_dim_2; ++n2) {
				const int d1 = (x1 * n1) % max_dim_1;
				const int d2 = (x2 * n2) % max_dim_2;
				root_of_unity[n1][n2] = Comp::GetPolar(
					RealScalar(1.0),
					RealScalar(2.0) * pi * (RealScalar(double(d1)) / RealScalar(double(max_dim_1)) + RealScalar(double(d2)) / RealScalar(double(max_dim_2))));
			}
		}

		const auto GetRoot = [&](const int i, const int j) {
			assert(0 <= i && i < max_dim_1 && 0 <= j && j < max_dim_2);
			return root_of_unity[i][j];
		};

		typedef std::pair<bool, Comp>MemComp;

		std::vector<MemComp>Z((n + 1) * 2, std::make_pair(false, Comp(0.0, 0.0)));
		std::vector<MemComp>Z1((n + 1) * (options.max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));
		std::vector<MemComp>Zb((n + 1) * (options.max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));
		std::vector<MemComp>Zm((n + 1) * (options.max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));
		std::vector<MemComp>Zm1((n + 1) * (options.max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));

		std::function<Comp(int, int)> GetZ;
		std::function<Comp(int, int)> GetZ1;
		std::function<Comp(int, int)> GetZb;
		std::function<Comp(int, int)> GetZm;
		std::function<Comp(int, int)> GetZm1;

		GetZ = [&](const int i, const int j) {
			assert(1 <= i && (i <= j || j == i - 1) && j <= n);
			if (i == j)return Comp(1.0, 0.0);
			if (j == i - 1)return Comp(1.0, 0.0);
			assert(i == 1 || j == n);
			const int index = (i == 1) ? j : (i + n + 1);
			assert(1 <= index && index < (n + 1) * 2);
			if (Z[index].first)return Z[index].second;
			Comp ans(0.0, 0.0);

			if (i == 1) {

				//no base pair
				ans += GetRoot(gz1_1(i, j), gz1_2(i, j));

				//there is only one outermost base pair (i, *)
				ans += GetZ1(i, j);

				//there are one or more outermost base pairs. The rightmost outermost base pair is (k + 1, *)
				for (int k = i; k <= j - 1; ++k) {
					const int r1 = gz0_1(i, j) - gz0_1(i, k) - gz0_1(k + 1, j);
					const int r2 = gz0_2(i, j) - gz0_2(i, k) - gz0_2(k + 1, j);
					ans += GetZ(i, k) * GetZ1(k + 1, j) * GetRoot(r1, r2);
				}
			}
			else {

				//i is not paired
				ans += GetZ(i + 1, j) * GetRoot(gz0_1(i, j) - gz0_1(i + 1, j), gz0_2(i, j) - gz0_2(i + 1, j));

				//(i,k) form a base pair
				for (int k = i + TURN + 1; k <= j && k - i <= options.max_span; ++k) {
					const int type = parasor_param::GetPairType(options.sequence[i - 1], options.sequence[k - 1]);
					if (type == 0)continue;
					const int r1 = gz0_1(i, j) - gz0_1(i, k) - gz0_1(k + 1, j);
					const int r2 = gz0_2(i, j) - gz0_2(i, k) - gz0_2(k + 1, j);
					ans += GetZb(i, k)
						* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (k - 1) + 1, true, options.sequence)))
						* GetZ(k + 1, j)
						* GetRoot(r1, r2);
				}
			}
			return (Z[index] = std::make_pair(true, ans)).second;
		};
		GetZ1 = [&](const int i, const int j) {
			assert(1 <= i && i <= j && j <= n);
			if (i == j)return Comp(0.0, 0.0);

			if ((j - i) > options.max_span) {
				return GetZ1(i, i + options.max_span) * GetRoot(gz3_1(i, j, i + options.max_span), gz3_2(i, j, i + options.max_span));
			}

			if (Z1[at(i, j)].first)return Z1[at(i, j)].second;

			Comp ans(0.0, 0.0);
			const int type = parasor_param::GetPairType(options.sequence[i - 1], options.sequence[j - 1]);
			if (type != 0) {
				ans += GetZb(i, j)
					* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, true, options.sequence)));
			}
			ans += GetZ1(i, j - 1) * GetRoot(gz3_1(i, j, j - 1), gz3_2(i, j, j - 1));

			return (Z1[at(i, j)] = std::make_pair(true, ans)).second;
		};
		GetZb = [&](const int i, const int j) {
			assert(1 <= i && i <= j && j <= n);
			if (i == j)return Comp(0.0, 0.0);

			const int type = parasor_param::GetPairType(options.sequence[i - 1], options.sequence[j - 1]);
			if (type == 0 || i + TURN >= j || (j - i) > options.max_span) {
				return Comp(0.0, 0.0);
			}

			if (Zb[at(i, j)].first)return Zb[at(i, j)].second;

			Comp ans(0.0, 0.0);

			//hairpin loop
			ans += Comp::LogRealToComp(RealScalar(parasor_param::ParHairpinEnergy(i - 1, j - 1, options.sequence))) * GetRoot(gz4_1(i, j), gz4_2(i, j));

			//internal loop
			for (int k = i + 1; k <= std::min(i + options.max_loop + 1, j - TURN - 2); ++k) {
				const int unpaired_base1 = k - i - 1;
				for (int l = std::max(k + TURN + 1, j - 1 - options.max_loop + unpaired_base1); l < j; ++l) {
					const int type = parasor_param::GetPairType(options.sequence[k - 1], options.sequence[l - 1]);
					if (type == 0)continue;
					ans += GetZb(k, l)
						* Comp::LogRealToComp(RealScalar(parasor_param::ParLoopEnergy(i - 1, j - 1, k - 1, l - 1, options.sequence)))
						* GetRoot(gz5_1(i, j, k, l), gz5_2(i, j, k, l));
				}
			}

			//multi loop
			const int rtype = parasor_param::GetPairTypeReverse(options.sequence[i - 1], options.sequence[j - 1]);
			for (int k = i + TURN + 3; k + TURN + 1 <= j - 1; ++k) {
				ans += GetZm(i + 1, k - 1)
					* GetZm1(k, j - 1)
					* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(rtype, (j - 1) - 1, (i - 1) + 1, false, options.sequence)))
					* Comp::LogRealToComp(RealScalar(parasor_param::ParMultiloopClosing()))
					* GetRoot(gz6_1(i, j, k), gz6_2(i, j, k));
			}

			return (Zb[at(i, j)] = std::make_pair(true, ans)).second;
		};
		GetZm = [&](const int i, const int j) {
			assert(1 <= i && (i <= j || j == i - 1) && j <= n);
			if (i == j)return Comp(0.0, 0.0);
			if (j == i - 1)return Comp(0.0, 0.0);
			if (Zm[at(i, j)].first)return Zm[at(i, j)].second;

			Comp ans(0.0, 0.0);
			for (int k = i; k + TURN + 1 <= j; ++k) {
				ans += (GetRoot(gz7_1(i, j, k), gz7_2(i, j, k)) + GetZm(i, k - 1) * GetRoot(gz8_1(i, j, k), gz8_2(i, j, k))) * GetZm1(k, j);
			}
			return (Zm[at(i, j)] = std::make_pair(true, ans)).second;
		};
		GetZm1 = [&](const int i, const int j) {
			assert(1 <= i && i <= j && j <= n);
			if (i == j)return Comp(0.0, 0.0);
			if (Zm1[at(i, j)].first)return Zm1[at(i, j)].second;

			Comp ans(0.0, 0.0);
			const int type = parasor_param::GetPairType(options.sequence[i - 1], options.sequence[j - 1]);
			if (type != 0) {
				ans += GetZb(i, j)
					* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, false, options.sequence)))
					* Comp::LogRealToComp(RealScalar(parasor_param::ParMultiloopInternal()));
			}
			ans += GetZm1(i, j - 1) * GetRoot(gz9_1(i, j, j - 1), gz9_2(i, j, j - 1));

			return (Zm1[at(i, j)] = std::make_pair(true, ans)).second;
		};

		//[Hagio et al., 2018], Algorithm2
		std::vector<MemComp>Wb((n + 1) * (options.max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));
		std::function<Comp(int, int)> GetWb;
		GetWb = [&](const int i, const int j) {
			assert(1 <= i && i <= j && j <= n);

			const int type = parasor_param::GetPairType(options.sequence[i - 1], options.sequence[j - 1]);
			if (type == 0 || i + TURN >= j || (j - i) > options.max_span) {
				return Comp(0.0, 0.0);
			}

			if (Wb[at(i, j)].first)return Wb[at(i, j)].second;

			Comp ans(0.0, 0.0);

			//(i,j) is an outermost base pair
			ans += GetZ(1, i - 1)
				* GetZ(j + 1, n)
				* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, true, options.sequence)))
				* GetRoot(gw1_1(i, j), gw1_2(i, j));

			//base pair (h,l) exists outside (i,j)
			for (int h = std::max(1, std::max(i - options.max_span + TURN + 2, i - options.max_loop - 1)); h < i; ++h) {
				for (int l = j + 1; l <= n && (l - h) <= options.max_span && (i - h - 1) + (l - j - 1) <= options.max_loop; ++l) {
					//internal loop, bulge, stem
					const int type = parasor_param::GetPairType(options.sequence[h - 1], options.sequence[l - 1]);
					if (type == 0)continue;
					ans += GetWb(h, l)
						* Comp::LogRealToComp(RealScalar(parasor_param::ParLoopEnergy(h - 1, l - 1, i - 1, j - 1, options.sequence)))
						* GetRoot(gw2_1(i, j, h, l), gw2_2(i, j, h, l));
				}
			}
			for (int h = std::max(1, i - options.max_span + TURN + 2); h < i; ++h) {
				for (int l = j + 1; l <= n && (l - h) <= options.max_span; ++l) {
					//(h,l) is multiloop closing base pair. (i,j) is one of the branches
					const int rtype = parasor_param::GetPairTypeReverse(options.sequence[h - 1], options.sequence[l - 1]);
					if (rtype == 0)continue;
					ans += Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, false, options.sequence)))
						* Comp::LogRealToComp(RealScalar(parasor_param::ParMultiloopInternal()))
						* GetWb(h, l)
						* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(rtype, (l - 1) - 1, (h - 1) + 1, false, options.sequence)))
						* Comp::LogRealToComp(RealScalar(parasor_param::ParMultiloopClosing()))
						* (
							GetZm(h + 1, i - 1) * GetRoot(gw3_1(i, j, h, l), gw3_2(i, j, h, l))
							+ GetZm(j + 1, l - 1) * GetRoot(gw4_1(i, j, h, l), gw4_2(i, j, h, l))
							+ GetZm(h + 1, i - 1) * GetZm(j + 1, l - 1) * GetRoot(gw5_1(i, j, h, l), gw5_2(i, j, h, l))
							);
				}
			}
			return (Wb[at(i, j)] = std::make_pair(true, ans)).second;
		};

		zeta[x1][x2] = GetZ(1, n);

		for (int i = 1; i <= n; ++i)for (int j = i; j <= n && (j - i) <= options.max_span; ++j) {

			//[Hagio et al., 2018].(9)
			prob[at(i, j)][x1][x2] = GetZb(i, j) * GetWb(i, j);
		}
	}

	std::vector<std::vector<std::vector<std::vector<Comp>>>>P(
		options.max_dim1 + 1, std::vector<std::vector<std::vector<Comp>>>(
			options.max_dim2 + 1, std::vector<std::vector<Comp>>(
				n + 1, std::vector<Comp>(
					options.max_span + 1, Comp(0.0, 0.0)))));

	std::vector<std::pair<int, int>>index;
	for (int i = 1; i <= n; ++i) for (int j = i + 1; j <= n && (j - i) <= options.max_span; ++j)index.push_back(std::make_pair(i, j));

	int count = 0;
#pragma omp parallel for schedule(dynamic, 1)
	for (int i = 0; i < int(index.size()); ++i) {
		const int a = index[i].first;
		const int b = index[i].second;

		std::vector<std::vector<Comp>>p_dash(max_dim_1, std::vector<Comp>(max_dim_2, Comp(0.0, 0.0)));
		for (int x1 = 0; x1 < max_dim_1; ++x1) {
			std::vector<Comp>p_temp_a(max_dim_2, Comp(0.0, 0.0));
			for (int x2 = 0; x2 < max_dim_2; ++x2)p_temp_a[x2] = prob[at(a, b)][x1][x2];
			const auto p_temp_b = FourierTransform(p_temp_a, options.allow_fft);
			for (int x2 = 0; x2 < max_dim_2; ++x2) {
				p_dash[x1][x2] = p_temp_b[x2];
			}
		}

		std::vector<std::vector<Comp>>p_(max_dim_1, std::vector<Comp>(max_dim_2, Comp(0.0, 0.0)));
		for (int x2 = 0; x2 < max_dim_2; ++x2) {
			std::vector<Comp>p_temp_a(max_dim_1, Comp(0.0, 0.0));
			for (int x1 = 0; x1 < max_dim_1; ++x1)p_temp_a[x1] = p_dash[x1][x2];
			const auto p_temp_b = FourierTransform(p_temp_a, options.allow_fft);
			for (int x1 = 0; x1 < max_dim_1; ++x1) {
				p_[x1][x2] = p_temp_b[x1];
			}
		}

		for (int d1 = 0; d1 <= options.max_dim1; ++d1)for (int d2 = 0; d2 <= options.max_dim2; ++d2) {
			P[d1][d2][a][b - a] = p_[d1][d2];
		}
	}

	std::vector<std::vector<Comp>>z_dash(max_dim_1, std::vector<Comp>(max_dim_2, Comp(0.0, 0.0)));
	for (int x1 = 0; x1 < max_dim_1; ++x1) {
		std::vector<Comp>z_temp_a(max_dim_2, Comp(0.0, 0.0));
		for (int x2 = 0; x2 < max_dim_2; ++x2)z_temp_a[x2] = zeta[x1][x2];
		const auto z_temp_b = FourierTransform(z_temp_a, options.allow_fft);
		for (int x2 = 0; x2 < max_dim_2; ++x2) {
			z_dash[x1][x2] = z_temp_b[x2];
		}
	}

	std::vector<std::vector<Comp>>z(max_dim_1, std::vector<Comp>(max_dim_2, Comp(0.0, 0.0)));
	for (int x2 = 0; x2 < max_dim_2; ++x2) {
		std::vector<Comp>z_temp_a(max_dim_1, Comp(0.0, 0.0));
		for (int x1 = 0; x1 < max_dim_1; ++x1)z_temp_a[x1] = z_dash[x1][x2];
		const auto z_temp_b = FourierTransform(z_temp_a, options.allow_fft);
		for (int x1 = 0; x1 < max_dim_1; ++x1) {
			z[x1][x2] = z_temp_b[x1];
		}
	}

	if (options.allow_fft && (max_dim_1 != options.max_dim1 || max_dim_2 != options.max_dim2)) {
		std::vector<std::vector<Comp>>zz(options.max_dim1 + 1, std::vector<Comp>(options.max_dim2 + 1, Comp(0.0, 0.0)));
		for (int d1 = 0; d1 <= options.max_dim1; ++d1)for (int d2 = 0; d2 <= options.max_dim2; ++d2)zz[d1][d2] = z[d1][d2];
		return std::make_pair(zz, P);
	}

	return std::make_pair(z, P);
}
std::vector<double> BruteForceRintD1Dim(const RintX1DOptions& options) {

	parasor_param::InitializeParameter(options.param_file_name, options.temperature);
	const int n = int(options.sequence.size());
	const std::vector<std::string> structures = EnumerateStructures(options.sequence, options.max_span, options.max_loop);

	std::vector<double>z(options.max_dim1 + 1, 0.0);

	for (const std::string s : structures) {
		const int d = ComputeHammingDistance(s, options.reference_structure1);
		const double e = EvalSpecificStructure(options.sequence, s);
		z[d] += e;
	}

	return z;
}
std::vector<double> BruteForceRintD1DimPK(const RintX1DOptions& options) {

	parasor_param::InitializeParameter(options.param_file_name, options.temperature);
	const int n = int(options.sequence.size());
	const std::vector<std::string> structures = EnumerateStructures(options.sequence, options.max_span, options.max_loop);

	std::vector<double>z(options.max_dim1 + 1, 0.0);

	for (const std::string s : structures) {
		const int d = ComputeHammingDistance(VerifyAndParseStructure(s, options.sequence, options.max_span, options.max_loop), options.S1);
		const double e = EvalSpecificStructure(options.sequence, s);
		z[d] += e;
	}

	return z;
}

std::vector<std::vector<double>> BruteForceRintD2Dim(const RintX2DOptions& options) {

	parasor_param::InitializeParameter(options.param_file_name, options.temperature);
	const int n = int(options.sequence.size());
	const std::vector<std::string> structures = EnumerateStructures(options.sequence, options.max_span, options.max_loop);

	std::vector<std::vector<double>>z(options.max_dim1 + 1, std::vector<double>(options.max_dim2 + 1, 0.0));

	for (const std::string s : structures) {
		const int d1 = ComputeHammingDistance(s, options.reference_structure1);
		const int d2 = ComputeHammingDistance(s, options.reference_structure2);
		const double e = EvalSpecificStructure(options.sequence, s);
		z[d1][d2] += e;
	}

	return z;
}
std::pair<std::vector<double>, std::vector<std::vector<std::vector<double>>>> BruteForceRintW1Dim(const RintX1DOptions& options) {

	parasor_param::InitializeParameter(options.param_file_name, options.temperature);
	const int n = int(options.sequence.size());
	const std::vector<std::string> structures = EnumerateStructures(options.sequence, options.max_span, options.max_loop);

	std::vector<double>z(options.max_dim1 + 1, 0);
	std::vector<std::vector<std::vector<double>>>ans(
		options.max_dim1 + 1, std::vector<std::vector<double>>(
			n + 1, std::vector<double>(
				options.max_span + 1, 0.0)));

	for (const std::string s : structures) {
		const int d = ComputeHammingDistance(s, options.reference_structure1);
		const double e = EvalSpecificStructure(options.sequence, s);
		z[d] += e;
		std::stack<int> stk;
		for (int i = 0; i < n; ++i) {
			switch (s[i]) {
			case '(':
				stk.push(i);
				break;
			case ')':
				ans[d][stk.top() + 1][i - stk.top()] += e;
				stk.pop();
				break;
			case '.':
				break;
			default:
				assert(0);
				break;
			}
		}
	}

	return std::make_pair(z, ans);
}
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<std::vector<std::vector<double>>>>> BruteForceRintW2Dim(const RintX2DOptions& options) {

	parasor_param::InitializeParameter(options.param_file_name, options.temperature);
	const int n = int(options.sequence.size());
	const std::vector<std::string> structures = EnumerateStructures(options.sequence, options.max_span, options.max_loop);

	std::vector<std::vector<double>>z(
		options.max_dim1 + 1, std::vector<double>(
			options.max_dim2 + 1, 0.0));

	std::vector<std::vector<std::vector<std::vector<double>>>>ans(
		options.max_dim1 + 1, std::vector<std::vector<std::vector<double>>>(
			options.max_dim2 + 1, std::vector<std::vector<double>>(
				n + 1, std::vector<double>(
					options.max_span + 1, 0.0))));

	for (const std::string s : structures) {
		const int d1 = ComputeHammingDistance(s, options.reference_structure1);
		const int d2 = ComputeHammingDistance(s, options.reference_structure2);
		const double e = EvalSpecificStructure(options.sequence, s);
		z[d1][d2] += e;
		std::stack<int> stk;
		for (int i = 0; i < n; ++i) {
			switch (s[i]) {
			case '(':
				stk.push(i);
				break;
			case ')':
				ans[d1][d2][stk.top() + 1][i - stk.top()] += e;
				stk.pop();
				break;
			case '.':
				break;
			default:
				assert(0);
				break;
			}
		}
	}

	return std::make_pair(z, ans);
}

std::pair<std::vector<Floating>, std::vector<std::vector<std::vector<Floating>>>>RintW1DimToBppm(
	const std::vector<WideComplexNumber<Floating>> z,
	const std::vector<std::vector<std::vector<WideComplexNumber<Floating>>>> p) {

	//ComputeHagioの返り値はlog scaleで、かつ"内側分配関数"と"塩基対確率行列の成分の分子"に分かれた状態なので、
	//これを塩基対確率行列の形に整える。

	//入力: ComputeHagio1dimの返り値の、firstがzでsecondがpとする。
	//返り値: firstはハミング距離ごとの存在確率で、secondはハミング距離ごとの塩基対確率行列。

	const int dim = int(z.size()) - 1;
	const int n = int(p[0].size()) - 1;
	const int max_span = int(p[0][0].size()) - 1;

	std::vector<Floating>z_ans_log(dim + 1, -std::numeric_limits<Floating>::infinity());
	Floating z_ans_log_max = -std::numeric_limits<Floating>::max();
	for (int x = 0; x <= dim; ++x) {
		if (z[x].real <= 0.0)continue;
		z_ans_log[x] = z[x].log_scale + log(z[x].real);
		z_ans_log_max = std::max(z_ans_log_max, z_ans_log[x]);
	}
	std::vector<Floating>z_ans(dim + 1, Floating(0.0));
	for (int x = 0; x <= dim; ++x)z_ans[x] = exp(z_ans_log[x] - z_ans_log_max);
	Floating z_sum = Floating(0.0);
	for (int x = 0; x <= dim; ++x)z_sum += z_ans[x];
	for (int x = 0; x <= dim; ++x)z_ans[x] /= z_sum;

	std::vector<std::vector<std::vector<Floating>>>p_ans(
		dim + 1, std::vector<std::vector<Floating>>(
			n + 1, std::vector<Floating>(
				max_span + 1, Floating(0.0))));
	for (int x = 0; x <= dim; ++x) {
		if (z_ans[x] == 0.0)continue;
		for (int i = 1; i <= n; ++i)for (int j = i + 1; j <= n && j - i <= max_span; ++j) {
			if (p[x][i][j - i].real <= 0.0)continue;
			const Floating p_log = p[x][i][j - i].log_scale + log(p[x][i][j - i].real);
			const Floating p_num = exp(p_log - z_ans_log[x]);
			p_ans[x][i][j - i] = std::min(Floating(1.0), p_num);
		}
	}

	return std::make_pair(z_ans, p_ans);
}

std::pair<std::vector<IntervalVar>, std::vector<std::vector<std::vector<IntervalVar>>>>RintW1DimToBppm(
	const std::vector<WideComplexNumber<IntervalVar>> z,
	const std::vector<std::vector<std::vector<WideComplexNumber<IntervalVar>>>> p) {

	//ComputeHagioの返り値はlog scaleで、かつ内側分配関数と塩基対確率行列の分子に分かれた状態なので、
	//これを塩基対確率行列の形に整える。

	const int dim = int(z.size()) - 1;
	const int n = int(p[0].size()) - 1;
	const int max_span = int(p[0][0].size()) - 1;
	const Floating FINF = std::numeric_limits<Floating>::infinity();

	std::vector<IntervalVar>z_ans_log(dim + 1, IntervalVar(-FINF, -FINF));
	Floating z_ans_log_max = -std::numeric_limits<Floating>::max();
	for (int x = 0; x <= dim; ++x) {
		assert(z[x].log_scale.upper() <= std::numeric_limits<Floating>::max());
		if (z[x].real.upper() <= 0.0)continue;
		const IntervalVar zz = /*kv::*/max(z[x].real, IntervalVar(0.0, 0.0));
		z_ans_log[x] = z[x].log_scale + log(zz);
		if (zz.lower() > 0.0)z_ans_log_max = std::max(z_ans_log_max, Floating(2.0 * z_ans_log[x].upper() - mid(z_ans_log[x])));
	}

	std::vector<IntervalVar>z_ans(dim + 1, IntervalVar(0.0, 0.0));
	for (int x = 0; x <= dim; ++x)z_ans[x] = exp(z_ans_log[x] - IntervalVar(z_ans_log_max));

	IntervalVar z_sum(0.0, 0.0);
	for (int x = 0; x <= dim; ++x)z_sum += z_ans[x];

	for (int x = 0; x <= dim; ++x) {
		z_ans[x] /= z_sum;
		z_ans[x] = /*kv::*/max(z_ans[x], IntervalVar(0.0, 0.0));
		z_ans[x] = /*kv::*/min(z_ans[x], IntervalVar(1.0, 1.0));
	}

	std::vector<std::vector<std::vector<IntervalVar>>>p_ans(
		dim + 1, std::vector<std::vector<IntervalVar>>(
			n + 1, std::vector<IntervalVar>(
				max_span + 1, IntervalVar(0.0, 0.0))));

	for (int x = 0; x <= dim; ++x) {
		if (z_ans[x].upper() == 0.0)continue;
		if (z_ans_log[x].upper() == -FINF)continue;
		for (int i = 1; i <= n; ++i)for (int j = i + 1; j <= n && j - i <= max_span; ++j) {
			assert(p[x][i][j - i].log_scale.upper() <= std::numeric_limits<Floating>::max());
			if (p[x][i][j - i].real.upper() <= 0.0)continue;
			const IntervalVar pp = /*kv::*/max(p[x][i][j - i].real, IntervalVar(0.0, 0.0));
			const IntervalVar p_log = p[x][i][j - i].log_scale + log(pp);
			if (p_log.upper() == -FINF)continue;

			p_ans[x][i][j - i] = exp(IntervalVar(p_log) - z_ans_log[x]);
			p_ans[x][i][j - i] = /*kv::*/max(p_ans[x][i][j - i], IntervalVar(0.0, 0.0));
			p_ans[x][i][j - i] = /*kv::*/min(p_ans[x][i][j - i], IntervalVar(1.0, 1.0));
		}
	}

	return std::make_pair(z_ans, p_ans);
}

std::vector<Floating>RegularizeRintD1Dim(const std::vector<WideComplexNumber<Floating>>& result) {

	//ComputeRintD1Dimの結果を実数にして、総和1に正規化して、logsumexpではない普通の数値型にして返す。

	typedef WideComplexNumber<Floating> Comp;
	const int dim1 = int(result.size());

	Comp sum(0.0, 0.0);
	for (int i = 0; i < dim1; ++i)if (result[i].real > Floating(0.0)) {
		sum += result[i];
	}

	//この時点で、結果の総和はexp(sum.log_scale)にほぼ等しい。

	std::vector<Floating> ans(dim1, Floating(0.0));
	for (int i = 0; i < dim1; ++i)if (result[i].real > Floating(0.0)) {
		Comp tmp;
		tmp.real = result[i].real;
		tmp.imag = result[i].imag;
		tmp.log_scale = result[i].log_scale - sum.log_scale;
		ans[i] = tmp.ToUsualComp().real;
	}

	Floating true_sum(0.0);
	for (int i = 0; i < dim1; ++i)if (result[i].real > Floating(0.0)) {
		assert(Floating(0.0) <= ans[i] && ans[i] <= Floating(1.0));
		true_sum += ans[i];
	}
	assert(Floating(0.0) < true_sum);
	for (int i = 0; i < dim1; ++i)if (result[i].real > Floating(0.0)) {
		ans[i] /= true_sum;
	}

	return ans;
}
std::vector<std::vector<Floating>>RegularizeRintD2Dim(const std::vector<std::vector<WideComplexNumber<Floating>>>& result) {

	//ComputeRintD2Dimの結果を実数にして、総和1に正規化して、logsumexpではない普通の数値型にして返す。

	typedef WideComplexNumber<Floating> Comp;
	const int dim1 = int(result.size());
	const int dim2 = int(result[0].size());

	Comp sum(0.0, 0.0);
	for (int i = 0; i < dim1; ++i)for (int j = 0; j < dim2; ++j)if (result[i][j].real > Floating(0.0)) {
		sum += result[i][j];
	}

	//この時点で、結果の総和はexp(sum.log_scale)にほぼ等しい。

	std::vector<std::vector<Floating>> ans(dim1, std::vector<Floating>(dim2, Floating(0.0)));
	for (int i = 0; i < dim1; ++i)for (int j = 0; j < dim2; ++j)if (result[i][j].real > Floating(0.0)) {
		Comp tmp;
		tmp.real = result[i][j].real;
		tmp.imag = result[i][j].imag;
		tmp.log_scale = result[i][j].log_scale - sum.log_scale;
		ans[i][j] = tmp.ToUsualComp().real;
	}

	Floating true_sum(0.0);
	for (int i = 0; i < dim1; ++i)for (int j = 0; j < dim2; ++j)if (result[i][j].real > Floating(0.0)) {
		assert(Floating(0.0) <= ans[i][j] && ans[i][j] <= Floating(1.0));
		true_sum += ans[i][j];
	}
	assert(Floating(0.0) < true_sum);
	for (int i = 0; i < dim1; ++i)for (int j = 0; j < dim2; ++j)if (result[i][j].real > Floating(0.0)) {
		ans[i][j] /= true_sum;
	}

	return ans;
}

std::pair<std::vector<std::vector<Floating>>, std::vector<std::vector<std::vector<std::vector<Floating>>>>>RintW2DimToBppm(
	const std::vector<std::vector<WideComplexNumber<Floating>>> z,
	const std::vector<std::vector<std::vector<std::vector<WideComplexNumber<Floating>>>>> p) {

	//ComputeHagioの返り値はlog scaleで、かつ"内側分配関数"と"塩基対確率行列の成分の分子"に分かれた状態なので、
	//これを塩基対確率行列の形に整える。

	//入力: ComputeHagio2dimの返り値の、firstがzでsecondがpとする。
	//返り値: firstはハミング距離ごとの存在確率で、secondはハミング距離ごとの塩基対確率行列。

	const int dim1 = int(z.size()) - 1;
	const int dim2 = int(z[0].size()) - 1;
	const int n = int(p[0][0].size()) - 1;
	const int max_span = int(p[0][0][0].size()) - 1;

	std::vector<std::vector<Floating>>z_ans_log(dim1 + 1, std::vector<Floating>(dim2 + 1, -std::numeric_limits<Floating>::infinity()));
	Floating z_ans_log_max = -std::numeric_limits<Floating>::max();
	for (int x1 = 0; x1 <= dim1; ++x1)for (int x2 = 0; x2 <= dim2; ++x2) {
		if (z[x1][x2].real <= 0.0)continue;
		z_ans_log[x1][x2] = z[x1][x2].log_scale + log(z[x1][x2].real);
		z_ans_log_max = std::max(z_ans_log_max, z_ans_log[x1][x2]);
	}
	std::vector<std::vector<Floating>>z_ans(dim1 + 1, std::vector<Floating>(dim2 + 1, Floating(0.0)));
	for (int x1 = 0; x1 <= dim1; ++x1)for (int x2 = 0; x2 <= dim2; ++x2)z_ans[x1][x2] = exp(z_ans_log[x1][x2] - z_ans_log_max);
	Floating z_sum = Floating(0.0);
	for (int x1 = 0; x1 <= dim1; ++x1)for (int x2 = 0; x2 <= dim2; ++x2)z_sum += z_ans[x1][x2];
	for (int x1 = 0; x1 <= dim1; ++x1)for (int x2 = 0; x2 <= dim2; ++x2)z_ans[x1][x2] /= z_sum;

	std::vector<std::vector<std::vector<std::vector<Floating>>>>p_ans(
		dim1 + 1, std::vector<std::vector<std::vector<Floating>>>(
			dim2 + 1, std::vector<std::vector<Floating>>(
				n + 1, std::vector<Floating>(
					max_span + 1, Floating(0.0)))));

	for (int x1 = 0; x1 <= dim1; ++x1)for (int x2 = 0; x2 <= dim2; ++x2) {
		if (z_ans[x1][x2] == 0.0)continue;
		for (int i = 1; i <= n; ++i)for (int j = i + 1; j <= n && j - i <= max_span; ++j) {
			if (p[x1][x2][i][j - i].real <= 0.0)continue;
			const Floating p_log = p[x1][x2][i][j - i].log_scale + log(p[x1][x2][i][j - i].real);
			const Floating p_num = exp(p_log - z_ans_log[x1][x2]);
			p_ans[x1][x2][i][j - i] = std::min(Floating(1.0), p_num);
		}
	}

	return std::make_pair(z_ans, p_ans);
}

int ComputeCredibilityLimit(const RintX1DOptions& options, const Floating threshold) {

	assert(Floating(0.0) <= threshold && threshold <= Floating(1.0));

	const int n = int(options.sequence.size());
	const auto tmp = ComputeRintD1Dim<WideComplexNumber<Floating>>(options);
	const auto mori_plot = RegularizeRintD1Dim(tmp);
	Floating sum = Floating(0.0);
	for (int i = 0; i <= n; ++i) {
		sum += mori_plot[i];
		if (sum >= threshold)return i;
	}

	return n;
}

std::vector<std::string> ComputeRintW1Dim_1(
	const RintX1DOptions& options,
	const int dim) {

	typedef WideComplexNumber<IntervalVar> Comp;
	typedef decltype(Comp::real) RealScalar;

	parasor_param::InitializeParameter(options.param_file_name, options.temperature);

	const int max_dim_ = options.allow_fft ? Ceiling2Power(options.max_dim1 + 1) : (options.max_dim1 + 1);

	const int n = int(options.sequence.size());
	const RealScalar pi = acos(RealScalar(-1.0));

	const std::vector<std::vector<int>>C = ComputePredistanceMatrix(options.S1);

	const auto gz1 = [&](const int i, const int j) {return C[i][j]; };
	const auto gz2 = [&](const int i, const int j, const int k) {return C[i][j] - C[i][k] - C[k + 1][j]; };
	const auto gz3 = [&](const int i, const int j, const int k) {return C[i][j] - C[i][k]; };
	const auto gz4 = [&](const int i, const int j) {return C[i][j] + 1 - 2 * options.S1[i][j]; };
	const auto gz5 = [&](const int i, const int j, const int k, const int l) {return C[i][j] - C[k][l] + 1 - 2 * options.S1[i][j]; };
	const auto gz6 = [&](const int i, const int j, const int k) {return C[i][j] - C[i + 1][k - 1] - C[k][j - 1] + 1 - 2 * options.S1[i][j]; };
	const auto gz7 = [&](const int i, const int j, const int k) {return C[i][j] - C[k][j]; };
	const auto gz8 = [&](const int i, const int j, const int k) {return C[i][j] - C[i][k - 1] - C[k][j]; };
	const auto gz9 = gz3;

	//[Hagio et al., 2018].(4)
	//Since there is a possibility that i == n + 1, an "if" statement is added.
	const auto gz0 = [&](const int i, const int j) {
		if (i > j)return 0;
		return C[i][j];
	};

	//[Hagio et al., 2018].(14)
	const auto gw1 = [&](const int i, const int j) {return gz0(1, n) - gz0(i, j) - gz0(1, i - 1) - gz0(j + 1, n); };
	const auto gw2 = [&](const int i, const int j, const int h, const int l) {return gz0(h, l) - gz0(i, j) + 1 - 2 * options.S1[h][l]; };
	const auto gw3 = [&](const int i, const int j, const int h, const int l) {return gw2(i, j, h, l) - gz0(h + 1, i - 1); };
	const auto gw4 = [&](const int i, const int j, const int h, const int l) {return gw2(i, j, h, l) - gz0(j + 1, l - 1); };
	const auto gw5 = [&](const int i, const int j, const int h, const int l) {return gw2(i, j, h, l) - gz0(h + 1, i - 1) - gz0(j + 1, l - 1); };

	const auto at = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n && j - i <= options.max_span);
		return i * (options.max_span + 1) + (j - i);
	};

	//std::vector<Comp>zeta(max_dim_, Comp(0.0, 0.0));
	//std::vector<std::vector<Comp>>prob((n + 1) * (options.max_span + 1), std::vector<Comp>(max_dim_, Comp(0.0, 0.0)));

	assert(0 <= dim && dim < max_dim_);
	const int x = dim;

	std::vector<Comp> root_of_unity(max_dim_);
	for (int n1 = 0; n1 < max_dim_; ++n1) {
		const int d = (x * n1) % max_dim_;
		root_of_unity[n1] = Comp::GetPolar(
			RealScalar(1.0),
			RealScalar(2.0) * pi * RealScalar(double(d)) / RealScalar(double(max_dim_)));
	}

	const auto GetRoot = [&](const int i) {
		assert(0 <= i && i < max_dim_);
		return root_of_unity[i];
	};

	typedef std::pair<bool, Comp>MemComp;

	std::vector<MemComp>Z((n + 1) * 2, std::make_pair(false, Comp(0.0, 0.0)));
	std::vector<MemComp>Z1((n + 1) * (options.max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));
	std::vector<MemComp>Zb((n + 1) * (options.max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));
	std::vector<MemComp>Zm((n + 1) * (options.max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));
	std::vector<MemComp>Zm1((n + 1) * (options.max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));

	std::function<Comp(int, int)> GetZ;
	std::function<Comp(int, int)> GetZ1;
	std::function<Comp(int, int)> GetZb;
	std::function<Comp(int, int)> GetZm;
	std::function<Comp(int, int)> GetZm1;

	GetZ = [&](const int i, const int j) {
		assert(1 <= i && (i <= j || j == i - 1) && j <= n);
		if (i == j)return Comp(1.0, 0.0);
		if (j == i - 1)return Comp(1.0, 0.0);
		assert(i == 1 || j == n);
		const int index = (i == 1) ? j : (i + n + 1);
		assert(1 <= index && index < (n + 1) * 2);
		if (Z[index].first)return Z[index].second;
		Comp ans(0.0, 0.0);

		if (i == 1) {

			//no base pair
			ans += GetRoot(gz1(i, j));

			//there is only one outermost base pair (i, *)
			ans += GetZ1(i, j);

			//there are one or more outermost base pairs. The rightmost outermost base pair is (k + 1, *)
			for (int k = i; k <= j - 1; ++k) {
				ans += GetZ(i, k) * GetZ1(k + 1, j) * GetRoot(/*gz2(i, j, k)*/gz0(i, j) - gz0(i, k) - gz0(k + 1, j));
			}
		}
		else {

			//i is not paired
			ans += GetZ(i + 1, j) * GetRoot(gz0(i, j) - gz0(i + 1, j));

			//(i,k) form a base pair
			for (int k = i + TURN + 1; k <= j && k - i <= options.max_span; ++k) {
				const int type = parasor_param::GetPairType(options.sequence[i - 1], options.sequence[k - 1]);
				if (type == 0)continue;
				ans += GetZb(i, k)
					* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (k - 1) + 1, true, options.sequence)))
					* GetZ(k + 1, j)
					* GetRoot(/*gz2(i, j, k)*/gz0(i, j) - gz0(i, k) - gz0(k + 1, j));
			}
		}
		return (Z[index] = std::make_pair(true, ans)).second;
	};
	GetZ1 = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n);
		if (i == j)return Comp(0.0, 0.0);

		if ((j - i) > options.max_span) {
			return GetZ1(i, i + options.max_span) * GetRoot(gz3(i, j, i + options.max_span));
		}

		if (Z1[at(i, j)].first)return Z1[at(i, j)].second;

		Comp ans(0.0, 0.0);
		const int type = parasor_param::GetPairType(options.sequence[i - 1], options.sequence[j - 1]);
		if (type != 0) {
			ans += GetZb(i, j)
				* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, true, options.sequence)));
		}
		ans += GetZ1(i, j - 1) * GetRoot(gz3(i, j, j - 1));

		return (Z1[at(i, j)] = std::make_pair(true, ans)).second;
	};
	GetZb = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n);
		if (i == j)return Comp(0.0, 0.0);

		const int type = parasor_param::GetPairType(options.sequence[i - 1], options.sequence[j - 1]);
		if (type == 0 || i + TURN >= j || (j - i) > options.max_span) {
			return Comp(0.0, 0.0);
		}

		if (Zb[at(i, j)].first)return Zb[at(i, j)].second;

		Comp ans(0.0, 0.0);

		//hairpin loop
		ans += Comp::LogRealToComp(RealScalar(parasor_param::ParHairpinEnergy(i - 1, j - 1, options.sequence))) * GetRoot(gz4(i, j));

		//internal loop
		for (int k = i + 1; k <= std::min(i + options.max_loop + 1, j - TURN - 2); ++k) {
			const int unpaired_base1 = k - i - 1;
			for (int l = std::max(k + TURN + 1, j - 1 - options.max_loop + unpaired_base1); l < j; ++l) {
				const int type = parasor_param::GetPairType(options.sequence[k - 1], options.sequence[l - 1]);
				if (type == 0)continue;
				ans += GetZb(k, l)
					* Comp::LogRealToComp(RealScalar(parasor_param::ParLoopEnergy(i - 1, j - 1, k - 1, l - 1, options.sequence)))
					* GetRoot(gz5(i, j, k, l));
			}
		}

		//multi loop
		const int rtype = parasor_param::GetPairTypeReverse(options.sequence[i - 1], options.sequence[j - 1]);
		for (int k = i + TURN + 3; k + TURN + 1 <= j - 1; ++k) {
			ans += GetZm(i + 1, k - 1)
				* GetZm1(k, j - 1)
				* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(rtype, (j - 1) - 1, (i - 1) + 1, false, options.sequence)))
				* Comp::LogRealToComp(RealScalar(parasor_param::ParMultiloopClosing()))
				* GetRoot(gz6(i, j, k));
		}

		return (Zb[at(i, j)] = std::make_pair(true, ans)).second;
	};
	GetZm = [&](const int i, const int j) {
		assert(1 <= i && (i <= j || j == i - 1) && j <= n);
		if (i == j)return Comp(0.0, 0.0);
		if (j == i - 1)return Comp(0.0, 0.0);
		if (Zm[at(i, j)].first)return Zm[at(i, j)].second;

		Comp ans(0.0, 0.0);
		for (int k = i; k + TURN + 1 <= j; ++k) {
			ans += (GetRoot(gz7(i, j, k)) + GetZm(i, k - 1) * GetRoot(gz8(i, j, k))) * GetZm1(k, j);
		}
		return (Zm[at(i, j)] = std::make_pair(true, ans)).second;
	};
	GetZm1 = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n);
		if (i == j)return Comp(0.0, 0.0);
		if (Zm1[at(i, j)].first)return Zm1[at(i, j)].second;

		Comp ans(0.0, 0.0);
		const int type = parasor_param::GetPairType(options.sequence[i - 1], options.sequence[j - 1]);
		if (type != 0) {
			ans += GetZb(i, j)
				* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, false, options.sequence)))
				* Comp::LogRealToComp(RealScalar(parasor_param::ParMultiloopInternal()));
		}
		ans += GetZm1(i, j - 1) * GetRoot(gz9(i, j, j - 1));

		return (Zm1[at(i, j)] = std::make_pair(true, ans)).second;
	};

	//[Hagio et al., 2018], Algorithm2
	std::vector<MemComp>Wb((n + 1) * (options.max_span + 1), std::make_pair(false, Comp(0.0, 0.0)));
	std::function<Comp(int, int)> GetWb;
	GetWb = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n);

		const int type = parasor_param::GetPairType(options.sequence[i - 1], options.sequence[j - 1]);
		if (type == 0 || i + TURN >= j || (j - i) > options.max_span) {
			return Comp(0.0, 0.0);
		}

		if (Wb[at(i, j)].first)return Wb[at(i, j)].second;

		Comp ans(0.0, 0.0);

		//(i,j) is an outermost base pair
		ans += GetZ(1, i - 1)
			* GetZ(j + 1, n)
			* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, true, options.sequence)))
			* GetRoot(gw1(i, j));

		//base pair (h,l) exists outside (i,j)
		for (int h = std::max(1, std::max(i - options.max_span + TURN + 2, i - options.max_loop - 1)); h < i; ++h) {
			for (int l = j + 1; l <= n && (l - h) <= options.max_span && (i - h - 1) + (l - j - 1) <= options.max_loop; ++l) {
				//internal loop, bulge, stem
				const int type = parasor_param::GetPairType(options.sequence[h - 1], options.sequence[l - 1]);
				if (type == 0)continue;
				ans += GetWb(h, l)
					* Comp::LogRealToComp(RealScalar(parasor_param::ParLoopEnergy(h - 1, l - 1, i - 1, j - 1, options.sequence)))
					* GetRoot(gw2(i, j, h, l));
			}
		}
		for (int h = std::max(1, i - options.max_span + TURN + 2); h < i; ++h) {
			for (int l = j + 1; l <= n && (l - h) <= options.max_span; ++l) {
				//(h,l) is multiloop closing base pair. (i,j) is one of the branches
				const int rtype = parasor_param::GetPairTypeReverse(options.sequence[h - 1], options.sequence[l - 1]);
				if (rtype == 0)continue;
				ans += Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(type, (i - 1) - 1, (j - 1) + 1, false, options.sequence)))
					* Comp::LogRealToComp(RealScalar(parasor_param::ParMultiloopInternal()))
					* GetWb(h, l)
					* Comp::LogRealToComp(RealScalar(parasor_param::ParDangling(rtype, (l - 1) - 1, (h - 1) + 1, false, options.sequence)))
					* Comp::LogRealToComp(RealScalar(parasor_param::ParMultiloopClosing()))
					* (
						GetZm(h + 1, i - 1) * GetRoot(gw3(i, j, h, l))
						+ GetZm(j + 1, l - 1) * GetRoot(gw4(i, j, h, l))
						+ GetZm(h + 1, i - 1) * GetZm(j + 1, l - 1) * GetRoot(gw5(i, j, h, l))
						);
			}
		}
		return (Wb[at(i, j)] = std::make_pair(true, ans)).second;
	};


	const auto dump = [](const Comp c) {
		std::string ans;
		union DU { double d; uint64_t u; };
		DU du;
		du.d = c.real.upper();
		ans += std::to_string(du.u) + std::string(",");
		du.d = c.real.lower();
		ans += std::to_string(du.u) + std::string(",");
		du.d = c.imag.upper();
		ans += std::to_string(du.u) + std::string(",");
		du.d = c.imag.lower();
		ans += std::to_string(du.u) + std::string(",");
		du.d = c.log_scale.upper();
		ans += std::to_string(du.u) + std::string(",");
		du.d = c.log_scale.lower();
		ans += std::to_string(du.u);
		return ans;
	};

	std::vector<std::string> ans;

	ans.push_back(dump(GetZ(1, n)));//zeta[x] = GetZ(1, n);

	for (int i = 1; i <= n; ++i)for (int j = i; j <= n && (j - i) <= options.max_span; ++j) {

		//[Hagio et al., 2018].(9)
		ans.push_back(dump(GetZb(i, j) * GetWb(i, j)));// prob[at(i, j)][x] = GetZb(i, j) * GetWb(i, j);
	}
	return ans;

//	//std::vector<std::vector<Comp>>prob((n + 1) * (options.max_span + 1), std::vector<Comp>(max_dim_, Comp(0.0, 0.0)));
//	std::vector<std::vector<std::vector<Comp>>>P(
//		options.max_dim1 + 1, std::vector<std::vector<Comp>>(
//			n + 1, std::vector<Comp>(
//				options.max_span + 1, Comp(0.0, 0.0))));
//
//	std::vector<std::pair<int, int>>index;
//	for (int i = 1; i <= n; ++i) for (int j = i + 1; j <= n && (j - i) <= options.max_span; ++j)index.push_back(std::make_pair(i, j));
//
//	int count = 0;
//#pragma omp parallel for schedule(dynamic, 1)
//	for (int i = 0; i < int(index.size()); ++i) {
//
//#pragma omp critical
//		{
//			if (i % 100 == 0) std::cout << i << "/" << int(index.size()) << std::endl;
//		}
//
//		const int a = index[i].first;
//		const int b = index[i].second;
//		const std::vector<Comp>p_ = FourierTransform(prob[at(a, b)], options.allow_fft);
//		for (int d = 0; d <= options.max_dim1; ++d)P[d][a][b - a] = p_[d];
//	}
//
//	const std::vector<Comp>z = FourierTransform(zeta, options.allow_fft);
//
//	if (options.allow_fft && max_dim_ != options.max_dim1) {
//		std::vector<Comp>zz(options.max_dim1 + 1, Comp(0.0, 0.0));
//		for (int d = 0; d <= options.max_dim1; ++d)zz[d] = z[d];
//		return std::make_pair(zz, P);
//	}
//
//	return std::make_pair(z, P);
}

std::vector<std::string> ComputeRintW1Dim_2(
	const RintX1DOptions& options,
	const std::vector<WideComplexNumber<IntervalVar>>& zeta,
	const std::vector<std::vector<WideComplexNumber<IntervalVar>>>& prob,
	const int num,
	const int mod) {

	typedef WideComplexNumber<IntervalVar> Comp;
	typedef decltype(Comp::real) RealScalar;

	parasor_param::InitializeParameter(options.param_file_name, options.temperature);

	const int max_dim_ = options.allow_fft ? Ceiling2Power(options.max_dim1 + 1) : (options.max_dim1 + 1);

	const int n = int(options.sequence.size());

	const auto at = [&](const int i, const int j) {
		assert(1 <= i && i <= j && j <= n && j - i <= options.max_span);
		return i * (options.max_span + 1) + (j - i);
	};

	const auto dump = [](const Comp c) {
		std::string ans;
		union DU { double d; uint64_t u; };
		DU du;
		du.d = c.real.upper();
		ans += std::to_string(du.u) + std::string(",");
		du.d = c.real.lower();
		ans += std::to_string(du.u) + std::string(",");
		du.d = c.imag.upper();
		ans += std::to_string(du.u) + std::string(",");
		du.d = c.imag.lower();
		ans += std::to_string(du.u) + std::string(",");
		du.d = c.log_scale.upper();
		ans += std::to_string(du.u) + std::string(",");
		du.d = c.log_scale.lower();
		ans += std::to_string(du.u);
		return ans;
	};

	std::vector<std::string> ans;
	//std::vector<std::vector<Comp>>prob((n + 1) * (options.max_span + 1), std::vector<Comp>(max_dim_, Comp(0.0, 0.0)));
	//std::vector<std::vector<std::vector<Comp>>>P(
	//	options.max_dim1 + 1, std::vector<std::vector<Comp>>(
	//		n + 1, std::vector<Comp>(
	//			options.max_span + 1, Comp(0.0, 0.0))));
	
	std::vector<std::pair<int, int>>index;
	for (int i = 1; i <= n; ++i) for (int j = i + 1; j <= n && (j - i) <= options.max_span; ++j)index.push_back(std::make_pair(i, j));

	for (int i = num; i < int(index.size()); i += mod) {
//	for (int i = 0; i < int(index.size()); ++i) {
		
		const int a = index[i].first;
		const int b = index[i].second;
		const std::vector<Comp>p_ = FourierTransform(prob[at(a, b)], options.allow_fft);
		for (int d = 0; d <= options.max_dim1; ++d) {
			//P[d][a][b - a] = p_[d];
			ans.push_back(dump(p_[d]));
		}
	}
	
	if (num == 0) {
		const std::vector<Comp>z = FourierTransform(zeta, options.allow_fft);

		if (options.allow_fft && max_dim_ != options.max_dim1) {
			//std::vector<Comp>zz(options.max_dim1 + 1, Comp(0.0, 0.0));
			for (int d = 0; d <= options.max_dim1; ++d)ans.push_back(dump(z[d]));//zz[d] = z[d];
		}
		else for (int d = 0; d < z.size(); ++d)ans.push_back(dump(z[d]));
	}
	return ans;

	//const std::vector<Comp>z = FourierTransform(zeta, options.allow_fft);
	//
	//if (options.allow_fft && max_dim_ != options.max_dim1) {
	//	std::vector<Comp>zz(options.max_dim1 + 1, Comp(0.0, 0.0));
	//	for (int d = 0; d <= options.max_dim1; ++d)zz[d] = z[d];
	//	return std::make_pair(zz, P);
	//}
	//
	//return std::make_pair(z, P);
}

}

