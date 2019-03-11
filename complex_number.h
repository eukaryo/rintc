/*
GNU GPL v2
Copyright (c) 2019 Hiroki Takizawa
*/

#ifndef RINTDWR_COMPLEX_NUMBER_H_
#define RINTDWR_COMPLEX_NUMBER_H_

#include <algorithm>
#include <cassert>

#include "interval_type.h"

namespace rintdwr {

template<typename RealScalar>class UsualComplexNumber {
public:
	RealScalar real, imag;
	UsualComplexNumber() :real(RealScalar(0.0)), imag(RealScalar(0.0)) {}
	UsualComplexNumber(const double r, const double i) :real(RealScalar(r)), imag(RealScalar(i)) {}

	UsualComplexNumber& operator += (const UsualComplexNumber& obj) {
		this->real = this->real + obj.real;
		this->imag = this->imag + obj.imag;
		return *this;
	}
	UsualComplexNumber& operator -= (const UsualComplexNumber& obj) {
		this->real = this->real - obj.real;
		this->imag = this->imag - obj.imag;
		return *this;
	}
	UsualComplexNumber& operator *= (const UsualComplexNumber& obj) {
		const RealScalar r = this->real * obj.real - this->imag * obj.imag;
		const RealScalar i = this->real * obj.imag + this->imag * obj.real;
		this->real = r;
		this->imag = i;
		return *this;
	}
	UsualComplexNumber& operator /= (const int& obj) {
		assert(obj > 0);
		RealScalar tmp(obj);
		this->real /= tmp;
		this->imag /= tmp;
		return *this;
	}
	UsualComplexNumber operator + (const UsualComplexNumber& obj)const { UsualComplexNumber re(*this); return re += obj; }
	UsualComplexNumber operator - (const UsualComplexNumber& obj)const { UsualComplexNumber re(*this); return re -= obj; }
	UsualComplexNumber operator * (const UsualComplexNumber& obj)const { UsualComplexNumber re(*this); return re *= obj; }
	UsualComplexNumber operator / (const int obj)const { UsualComplexNumber re(*this); return re /= obj; }

	static UsualComplexNumber GetPolar(const RealScalar& phase, const RealScalar& arg_radian) {
		UsualComplexNumber ans;
		ans.real = phase*cos(arg_radian);
		ans.imag = phase*sin(arg_radian);
		return ans;
	}
	UsualComplexNumber<RealScalar>ToUsualComp()const {
		UsualComplexNumber<RealScalar>ans;
		ans.real = this->real;
		ans.imag = this->imag;
		return ans;
	}
	static UsualComplexNumber LogRealToComp(const RealScalar& x) {
		//é¿êîxÇéÛÇØéÊÇ¡Çƒexp(x)+0iÇï‘Ç∑ÅB
		UsualComplexNumber ans;
		ans.real = exp(x);
		ans.imag = RealScalar(0.0);
		return ans;
	}
};

template<typename RealScalar>class WideComplexNumber {
private:
	static std::pair<IntervalVar, std::pair<IntervalVar, IntervalVar>> regularize_(const IntervalVar&, const IntervalVar&, const IntervalVar&);
	static std::pair<Floating, std::pair<Floating, Floating>> regularize_(const Floating&, const Floating&, const Floating&);

	static WideComplexNumber add_(
		const IntervalVar&,
		const IntervalVar&,
		const IntervalVar&,
		const IntervalVar&,
		const IntervalVar&,
		const IntervalVar&);
	static WideComplexNumber add_(
		const Floating&,
		const Floating&,
		const Floating&,
		const Floating&,
		const Floating&,
		const Floating&);

	static bool IsZero(const IntervalVar& re, const IntervalVar& im) {
		return (re * re + im * im).upper() == 0.0;
	}
	static bool IsZero(const Floating& re, const Floating& im) {
		return re * re + im * im == 0.0;
	}

public:
	static WideComplexNumber regularize(const RealScalar log_sc, const RealScalar re, const RealScalar im) {
		const auto s = regularize_(log_sc, re, im);
		WideComplexNumber ans;
		ans.log_scale = s.first;
		ans.real = s.second.first;
		ans.imag = s.second.second;
		return ans;
	}

	RealScalar log_scale;
	RealScalar real;
	RealScalar imag;

	WideComplexNumber() :log_scale(RealScalar(0.0)), real(RealScalar(0.0)), imag(RealScalar(0.0)) {}
	WideComplexNumber(const double r, const double i) {
		const auto a = regularize(RealScalar(0.0), RealScalar(r), RealScalar(i));

		this->log_scale = a.log_scale;
		this->real = a.real;
		this->imag = a.imag;
	}

	WideComplexNumber& operator += (const WideComplexNumber& obj) {

		if (IsZero(this->real, this->imag)) {
			this->log_scale = obj.log_scale;
			this->real = obj.real;
			this->imag = obj.imag;
			return *this;
		}
		if (IsZero(obj.real, obj.imag)) {
			return *this;
		}

		const auto p = add_(this->log_scale, this->real, this->imag, obj.log_scale, obj.real, obj.imag);
		this->log_scale = p.log_scale;
		this->real = p.real;
		this->imag = p.imag;

		return *this;
	}
	WideComplexNumber& operator -= (const WideComplexNumber& obj) {

		if (IsZero(this->real, this->imag)) {
			this->log_scale = obj.log_scale;
			this->real = -obj.real;
			this->imag = -obj.imag;
			return *this;
		}
		if (IsZero(obj.real, obj.imag)) {
			return *this;
		}

		const auto p = add_(this->log_scale, this->real, this->imag, obj.log_scale, -obj.real, -obj.imag);
		this->log_scale = p.log_scale;
		this->real = p.real;
		this->imag = p.imag;

		return *this;
	}
	WideComplexNumber& operator *= (const WideComplexNumber& obj) {

		if (IsZero(this->real, this->imag)) {
			return *this;
		}
		if (IsZero(obj.real, obj.imag)) {
			this->log_scale = RealScalar(0.0);
			this->real = RealScalar(0.0);
			this->imag = RealScalar(0.0);
			return *this;
		}

		const RealScalar r = (this->real * obj.real) - (this->imag * obj.imag);
		const RealScalar i = (this->real * obj.imag) + (this->imag * obj.real);
		this->log_scale += obj.log_scale;
		this->real = r;
		this->imag = i;
		return *this;
	}
	WideComplexNumber& operator /= (const int& obj) {

		assert(obj > 0);
		if (IsZero(this->real, this->imag)) {
			return *this;
		}

		this->log_scale -= log(RealScalar(obj));
		return *this;
	}
	WideComplexNumber operator + (const WideComplexNumber& obj)const { WideComplexNumber re(*this); return re += obj; }
	WideComplexNumber operator - (const WideComplexNumber& obj)const { WideComplexNumber re(*this); return re -= obj; }
	WideComplexNumber operator * (const WideComplexNumber& obj)const { WideComplexNumber re(*this); return re *= obj; }
	WideComplexNumber operator / (const int obj)const { WideComplexNumber re(*this); return re /= obj; }
	static WideComplexNumber GetPolar(const RealScalar& phase, const RealScalar& arg_radian) {
		WideComplexNumber ans;
		ans.log_scale = log(phase);
		ans.real = cos(arg_radian);
		ans.imag = sin(arg_radian);
		return ans;
	}
	bool IsZero()const {
		return WideComplexNumber::IsZero(this->real, this->imag);
	}
	UsualComplexNumber<RealScalar>ToUsualComp()const {
		UsualComplexNumber<RealScalar>ans;
		ans.real = this->real * exp(this->log_scale);
		ans.imag = this->imag * exp(this->log_scale);
		return ans;
	}
	static WideComplexNumber LogRealToComp(const RealScalar& x) {
		WideComplexNumber ans;
		ans.log_scale = RealScalar(x);
		ans.real = RealScalar(1.0);
		ans.imag = RealScalar(0.0);
		return ans;
	}
};

}

#endif//RINTDWR_COMPLEX_NUMBER_H_