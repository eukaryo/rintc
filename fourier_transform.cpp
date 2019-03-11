/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#include"fourier_transform.h"

#include"complex_number.h"
#include"interval_type.h"
#include"misc.h"

namespace rintdwr {

template<typename Comp>static std::vector<Comp>Dft(const std::vector<Comp>& s) {

	typedef decltype(Comp::real) RealScalar;

	const int n = int(s.size());
	const RealScalar pi = acos(RealScalar(-1.0));

	std::vector<Comp>ans(n);
	for (int i = 0; i < n; ++i) {
		Comp sum(0.0, 0.0);
		for (int j = 0; j < n; ++j) {
			const Comp freq = Comp::GetPolar(RealScalar(1.0),
				RealScalar(-2.0) * pi * RealScalar(double((i * j) % n)) / RealScalar(double(n)));
			sum += s[j] * freq;
		}
		ans[i] = sum / n;
	}

	return ans;
}
template<typename Comp>static std::vector<Comp>Fft(const std::vector<Comp>& s) {

	const auto bit_reverse = [](int x, int n) {
		int rev = 0;
		for (--n; n > 0; n >>= 1, x >>= 1)
			rev = (rev << 1) | (x & 1);
		return rev;
	};

	typedef decltype(Comp::real) RealScalar;

	const int N = int(s.size());
	assert(1 <= N && N <= 1000000);
	assert((N & (N - 1)) == 0);//N must be a power of 2.
	const RealScalar pi = acos(RealScalar(-1.0));
	std::vector<Comp> ans(N);
	for (int i = 0; i < N; ++i) {
		const int rev = bit_reverse(i, N);
		ans[i] = s[rev];
	}
	for (int b = 2; b <= N; b <<= 1) {
		for (int i = 0; i < N; ++i) {
			if (i % b >= b / 2)continue;
			const int k = i + b / 2;
			const RealScalar fi = RealScalar(-2.0) * pi * RealScalar(i % b) / RealScalar(b);
			const RealScalar fk = RealScalar(-2.0) * pi * RealScalar(k % b) / RealScalar(b);
			Comp Wi, Wk;
			Wi.real = cos(fi);
			Wi.imag = sin(fi);
			Wk.real = cos(fk);
			Wk.imag = sin(fk);
			ans[i] += Wi*ans[k];
			ans[k] = Wk*ans[k] + (ans[i] - Wi*ans[k]);
		}
	}
	for (int i = 0; i < N; ++i)ans[i] /= N;
	return ans;
}

std::vector<WideComplexNumber<IntervalVar>>FourierTransform(const std::vector<WideComplexNumber<IntervalVar>>& s, const bool allow_fft) {
	const int n = int(s.size());
	const bool can_use_fft = allow_fft && (n == Ceiling2Power(n));
	return can_use_fft ?
		Fft<WideComplexNumber<IntervalVar>>(s) :
		Dft<WideComplexNumber<IntervalVar>>(s);
}
std::vector<WideComplexNumber<Floating>>FourierTransform(const std::vector<WideComplexNumber<Floating>>& s, const bool allow_fft) {
	const int n = int(s.size());
	const bool can_use_fft = allow_fft && (n == Ceiling2Power(n));
	return can_use_fft ?
		Fft<WideComplexNumber<Floating>>(s) :
		Dft<WideComplexNumber<Floating>>(s);
}
std::vector<UsualComplexNumber<IntervalVar>>FourierTransform(const std::vector<UsualComplexNumber<IntervalVar>>& s, const bool allow_fft) {
	const int n = int(s.size());
	const bool can_use_fft = allow_fft && (n == Ceiling2Power(n));
	return can_use_fft ?
		Fft<UsualComplexNumber<IntervalVar>>(s) :
		Dft<UsualComplexNumber<IntervalVar>>(s);
}
std::vector<UsualComplexNumber<Floating>>FourierTransform(const std::vector<UsualComplexNumber<Floating>>& s, const bool allow_fft) {
	const int n = int(s.size());
	const bool can_use_fft = allow_fft && (n == Ceiling2Power(n));
	return can_use_fft ?
		Fft<UsualComplexNumber<Floating>>(s) :
		Dft<UsualComplexNumber<Floating>>(s);
}

}