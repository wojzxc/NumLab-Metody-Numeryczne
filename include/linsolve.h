#pragma once
#include <vector>
#include <stdexcept>   

namespace numlab {

	using Vector = std::vector<double>;
	using Matrix = std::vector<Vector>;

	Vector gaussian_elimination(Matrix A, Vector b, bool verbose = false);

} 
