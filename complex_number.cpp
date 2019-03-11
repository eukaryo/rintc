/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#include"complex_number.h"

#include <algorithm>
#include <cassert>

#include"interval_type.h"

namespace rintdwr {

//explicit instantiation

template std::pair<IntervalVar, std::pair<IntervalVar, IntervalVar>> WideComplexNumber<IntervalVar>::regularize_(const IntervalVar& log_sc, const IntervalVar& re, const IntervalVar& im);
template std::pair<Floating, std::pair<Floating, Floating>> WideComplexNumber<Floating>::regularize_(const Floating& log_sc, const Floating& re, const Floating& im);
template WideComplexNumber<IntervalVar> WideComplexNumber<IntervalVar>::add_(
	const IntervalVar& x_log_scale,
	const IntervalVar& x_real,
	const IntervalVar& x_imag,
	const IntervalVar& y_log_scale,
	const IntervalVar& y_real,
	const IntervalVar& y_imag);
template WideComplexNumber<Floating> WideComplexNumber<Floating>::add_(
	const Floating& x_log_scale,
	const Floating& x_real,
	const Floating& x_imag,
	const Floating& y_log_scale,
	const Floating& y_real,
	const Floating& y_imag);

//function definition

template<typename RealScalar>
std::pair<IntervalVar, std::pair<IntervalVar, IntervalVar>> WideComplexNumber<RealScalar>::regularize_(const IntervalVar& log_sc, const IntervalVar& re, const IntervalVar& im) {
	const Floating dist = (re * re + im * im).upper();
	if (dist == 0.0)return std::make_pair(IntervalVar(0.0), std::make_pair(IntervalVar(0.0), IntervalVar(0.0)));
	const IntervalVar t(1.0 / (sqrt(dist)));
	return std::make_pair(log_sc - log(t), std::make_pair(re*t, im*t));
}

template<typename RealScalar>
std::pair<Floating, std::pair<Floating, Floating>> WideComplexNumber<RealScalar>::regularize_(const Floating& log_sc, const Floating& re, const Floating& im) {
	const Floating dist = re * re + im * im;
	if (dist == 0.0)return std::make_pair(Floating(0.0), std::make_pair(Floating(0.0), Floating(0.0)));
	const Floating t(1.0 / (sqrt(dist)));
	return std::make_pair(log_sc - log(t), std::make_pair(re*t, im*t));
}

template<typename RealScalar>
WideComplexNumber<RealScalar> WideComplexNumber<RealScalar>::add_(
	const IntervalVar& x_log_scale,
	const IntervalVar& x_real,
	const IntervalVar& x_imag,
	const IntervalVar& y_log_scale,
	const IntervalVar& y_real,
	const IntervalVar& y_imag) {

	if (x_log_scale.upper() >= y_log_scale.upper()) {
		const IntervalVar p(x_log_scale.upper() - mid(x_log_scale) + y_log_scale.upper() - mid(y_log_scale));
		const IntervalVar c1 = exp(-p), c2 = exp(y_log_scale - x_log_scale - p);
		const IntervalVar re_ = c1 * x_real + c2 * y_real;
		const IntervalVar im_ = c1 * x_imag + c2 * y_imag;
		return regularize(x_log_scale + p, re_, im_);
	}
	const IntervalVar p(y_log_scale.upper() - mid(y_log_scale) + x_log_scale.upper() - mid(x_log_scale));
	const IntervalVar c1 = exp(-p), c2 = exp(x_log_scale - y_log_scale - p);
	const IntervalVar re_ = c1 * y_real + c2 * x_real;
	const IntervalVar im_ = c1 * y_imag + c2 * x_imag;
	return regularize(y_log_scale + p, re_, im_);
}

template<typename RealScalar>
WideComplexNumber<RealScalar> WideComplexNumber<RealScalar>::add_(
	const Floating& x_log_scale,
	const Floating& x_real,
	const Floating& x_imag,
	const Floating& y_log_scale,
	const Floating& y_real,
	const Floating& y_imag) {

	if (x_log_scale >= y_log_scale) {
		const Floating re_ = x_real + exp(y_log_scale - x_log_scale) * y_real;
		const Floating im_ = x_imag + exp(y_log_scale - x_log_scale) * y_imag;
		return regularize(x_log_scale, re_, im_);
	}
	const Floating re_ = y_real + exp(x_log_scale - y_log_scale) * x_real;
	const Floating im_ = y_imag + exp(x_log_scale - y_log_scale) * x_imag;
	return regularize(y_log_scale, re_, im_);
}

}