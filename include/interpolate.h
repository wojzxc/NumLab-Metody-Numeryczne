#pragma once
#include <vector>

namespace numlab {

	using Vector = std::vector<double>;
	using Matrix = std::vector<Vector>;

	double lagrange(const Vector& xi, const Vector& fi, double x);

	Vector newton_coeff(const Vector& xi, const Vector& fi);

	double newton_eval(const Vector& a, const Vector& xi, double x);

	double poly_eval(const Vector& a, double x);

	double poly_horner(const Vector& a, double x);

} 

