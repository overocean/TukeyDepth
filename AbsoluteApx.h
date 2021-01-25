#ifndef ABSOLUTEAPX_H_
#define ABSOLUTEAPX_H_

#include <limits>
#include <iostream>
#include <functional>
#include <random>
#include <unordered_set>
#include <thread>

#include <gmpxx.h>
#include <eigen3/Eigen/Dense>

#include "types.h"

namespace TD {

#define SAMPLE_LIMIT 5

/*
 * A Monte Carlo approximation algorithm for Tukey depth in arbitrary dimension.
 * This approximation has an absolute performance guarantee of bounded error
 * with probability no less than (1 - 1/e).
 * @warning the points in the data set need to be in general position.
 * @tparam DataType type for coordinates of the points in DataInfo
 * @tparam PrecisionType real number type for precision of calculation.
 */
template<typename DataType, typename PrecisionType>
class AbsoluteApx {

typedef std::function<SizeType(const PointSet<PrecisionType>&, const SizeType)> KDSolverType;

public:
	/* Constructor
	 * @param number_of_threads optional parameter for number of threads to use.
	 */
	explicit AbsoluteApx(const int number_of_threads = 0) {
		int thread_num = number_of_threads;
		if (thread_num <= 0) {
			thread_num = std::thread::hardware_concurrency();
			if (thread_num <= 0) {
				thread_num = 1;
			}
		}

		std::random_device rd;
		m_rngs.resize(thread_num);
		for (auto &rng: m_rngs) {
			rng.seed(rd());
		}

		m_results.resize(thread_num);
		m_solvers.resize(thread_num);
	}
	virtual ~AbsoluteApx(){};

	/* Function to do the approximation for a point.
	 * @param data_info the point set
	 * @param index the index of the point in the point set to compute depth for.
	 * @param sigma the bounded error
	 * @param k the dimension to reduce the problem to.
	 * @param solver the function to solve the k-dimensional problem. Ignored when k == 1,
	 * @param early_stop optional parameter to indicate whether we should stop the search
	 * after finding a result no more than sigma. This can shorten the running time for
	 * points with depth no more than sigma.
	 * @return the approximated depth.
	 */
	SizeType approxDepth(const PointSet<DataType>& data_info, const SizeType index,
			const SizeType sigma, const int k, KDSolverType solver = nullptr, const bool early_stop = false) {
		if (k >= data_info.m_dim || k <= 0) {
			throw std::runtime_error("Incorrect value for parameter k.");
		}

		if (k > 2 && solver == nullptr) {
			throw std::runtime_error("A solver for k-dimensional problem is required.");
		}

		if (data_info.m_point_num < data_info.m_dim + k) {
			std::cout << "This problem is too trivial. All points would have depth 0, if they are in general position.\n";
			return 0;
		}

		if (index >= data_info.m_point_num) {
			throw std::runtime_error("Bad point index!");
		}

		if (sigma >= data_info.m_point_num) {
			throw std::runtime_error("Bad sigma value!");
		}

		m_data = data_info.m_data;
		m_point_num = data_info.m_point_num;
		m_dim = data_info.m_dim;

		for (auto& thread_solver: m_solvers) {
			thread_solver = solver;
		}

		int num_threads = m_results.size();
		std::vector<std::thread> threads(num_threads);

		mpz_class sampleNumMpz = computeSampleNum(sigma, k);
		mpz_class ull_max(std::to_string(std::numeric_limits<unsigned long long int>::max()).c_str());

		for (int i = 0; i < num_threads; ++i) {
			mpz_class thread_sample_num = sampleNumMpz / num_threads;
			if (num_threads > 1 && i + 1 == num_threads) {
				thread_sample_num = sampleNumMpz - thread_sample_num * i;
			}

			if (thread_sample_num > ull_max) {
				threads[i] = std::thread( [this, i, index, thread_sample_num, sigma, k, early_stop] {
					approxDepthImpl(i, index, thread_sample_num, sigma, k, early_stop);
				});
			} else if (thread_sample_num > std::numeric_limits<unsigned int>::max()) {
				unsigned long long int sampleNum = 0;
				mpz_export(&sampleNum, 0, -1, sizeof sampleNum, 0, 0, thread_sample_num.get_mpz_t());
				threads[i] = std::thread( [this, i, index, sampleNum, sigma, k, early_stop] {
					approxDepthImpl(i, index, sampleNum, sigma, k, early_stop);
				});
			} else {
				unsigned int sampleNum = (unsigned int)thread_sample_num.get_ui();
				threads[i] = std::thread( [this, i, index, sampleNum, sigma, k, early_stop] {
					approxDepthImpl(i, index, sampleNum, sigma, k, early_stop);
				});
			}
		}

		for (int i = 0; i < num_threads; ++i) {
			threads[i].join();
		}

		SizeType depth = *std::min_element(m_results.begin(), m_results.end());

		return depth;
	}

private:
	DataType * m_data = 0;  // points stored as two dimensional array
	SizeType m_point_num = 0;  // number of points
	SizeType m_dim = 0;  // dimension

	// objects not shared by threads
	std::vector<std::mt19937> m_rngs; //random number generator
	std::vector<SizeType> m_results;
	std::vector<KDSolverType> m_solvers;

	template<typename SampleNumType>
	void approxDepthImpl(const int thread_id, const SizeType index, SampleNumType sampleNum,
			const SizeType sigma, const int k, const bool early_stop) {
		SizeType depth_min = m_point_num;  // the final depth
		std::uniform_int_distribution<> distrib(0, depth_min - 1);
		for (; sampleNum > 0; --sampleNum) {
			SizeType cur_depth = randomUpperBound(thread_id, index, k, distrib);
			depth_min = std::min(depth_min, cur_depth);

			if (depth_min == 0) break;  // depth can not be smaller
			if (early_stop && depth_min <= sigma) break;  // error of the depth is small enough
		}

		m_results[thread_id] = depth_min;
	}

	SizeType randomUpperBound(const int thread_id, const SizeType p, const int k,
			std::uniform_int_distribution<>& distrib) {

		//Sampling d-k points, we will have a (d-k)-flat defined by those d-k points and p.
		//Then we can find k-flat orthogonal to this (d-k)-flat, and project the point set S
		//on to the k-flat to get a k-dimensional problem, which can give an upper bound.

		//In this function we are going to find k orthogonal vectors in the k-flat, and compute
		//each point's coordinates by its dot product with the k vectors.
		std::vector<Eigen::Matrix<PrecisionType, Eigen::Dynamic, 1>> vectors(k);
		for (auto v: vectors) {
			v.resize(m_dim);
		}

		//We first take d-1 pints from the S, and find the k normal vectors of the (d-k)-flat
		//defined by d-k of them plus p.
		std::unordered_set<SizeType> sample_point_set;
		int sample_counter = 0;
		do {
			if (sample_counter > SAMPLE_LIMIT) {
				std::cerr << "Can't get samples points in general position after "
						  << SAMPLE_LIMIT << " attempts.\n";
				throw std::runtime_error("Data set is too degenerated.");
			}
			sample_point_set.clear();
			sample_point_set.insert(p);
			randomSample(thread_id, distrib, sample_point_set);
			++sample_counter;
		} while (!getNormalVectors(p, sample_point_set, vectors));

		SizeType projected_point_num = m_point_num - m_dim + k;
		std::vector<PrecisionType> data_set(projected_point_num * k);

		//mapping p to the k-flat
		for (int i = 0; i < k; ++i) {
			data_set[i] = innerProduct(vectors[i], m_data + p * m_dim);
		}

		//mapping S (exclude sample_point_set) to the k-flat
		SizeType idx = k;
		for (SizeType i = 0; i < m_point_num; ++i) {
			if (sample_point_set.count(i) == 0) {
				for (int j = 0; j < k; ++j) {
					data_set[idx] = innerProduct(vectors[j], m_data + i * m_dim);
					++idx;
				}
			}
		}

		// The last d-k points should not be counted for the depth value, since we can perturb
		// the halfspace a little bit to exclude those points from the halfspace.
		auto itr = sample_point_set.begin();
		for (int i = 1; i < k; ++i) {
			SizeType cur_point = *itr;
			for (int j = 0; j < k; ++j) {
				data_set[idx] = innerProduct(vectors[j], m_data + cur_point * m_dim);
				++idx;
			}
			++itr;
		}

		SizeType depth = 0;
		if (k == 1) {
			SizeType depth1 = 0;
			SizeType depth2 = 0;  // depth1 and depth2 are the number of points in the two closed halfspaces
			PrecisionType p_val = data_set[0];
			for (SizeType i = 1; i < projected_point_num; ++i) {
				if (data_set[i] <= p_val)  //depth is defined by a closed halfspace
					depth1++;
				if (data_set[i] >= p_val)
					depth2++;
			}
			depth = std::min(depth1, depth2);
		} else {
			PointSet<PrecisionType> kdim_set(data_set.data(), projected_point_num, k);
			depth = m_solvers[thread_id](kdim_set, 0);
		}

		return depth;
	}

	void randomSample(const int thread_id, std::uniform_int_distribution<>& distrib, std::unordered_set<SizeType> &result) {
		while (result.size() < m_dim) {
			result.insert(distrib(m_rngs[thread_id]));
		}
	}

	PrecisionType innerProduct(Eigen::Matrix<PrecisionType, Eigen::Dynamic, 1> &a, const DataType *b) {
		PrecisionType prod(0);
		for (std::size_t i = 0; i < m_dim; ++i ) {
				prod += a(i, 0) * b[ i ];
		}

		return prod;
	}

	/*
	 * Function to find k vectors that are orthogonal to the (d-k) flat defined by p and last d-k points
	 * in sample_point_set. The result k vectors are orthogonal to each other too.
	 */
	bool getNormalVectors(const SizeType p, const std::unordered_set<SizeType>& sample_point_set,
		std::vector<Eigen::Matrix<PrecisionType, Eigen::Dynamic, 1>>& result) const {

		Eigen::Matrix<PrecisionType, Eigen::Dynamic, Eigen::Dynamic> mat;
		mat.resize(m_dim, m_dim);
		Eigen::Matrix<PrecisionType, Eigen::Dynamic, 1> rhs;
		rhs.resize(m_dim);

		//The first equation is like 'x1 + x2 + ... + xd = d'
		for (SizeType i = 0; i < m_dim; ++i) {
			mat(0, i) = 1;
			rhs(i, 0) = 0;
		}
		rhs(0, 0) = m_dim;  // can cause overflow?

		//equation i is like '(a_i^1 -a_p^1)*x1 + (a_i^2 -a_p^2)*x2 + ... + (a_i^d -a_p^d)*xd = 0'
		SizeType i = 1;
		for (auto itr = sample_point_set.begin(); itr != sample_point_set.end(); ++itr) {
			SizeType cur_point = *itr;
			if (cur_point != p) {
				for (SizeType j = 0; j < m_dim; ++j) {
					mat(i, j) = m_data[cur_point * m_dim + j] - m_data[p * m_dim + j];
				}
				++i;
			}
		}

		SizeType vector_num = result.size();
		for (SizeType i = 0; i < vector_num; ++i) {
			//if (std::is_floating_point<PrecisionType>::value) {
			//	result[i] = mat.fullPivLu().solve(rhs);
			//} else {
			//	result[i] = mat.householderQr().solve(rhs);  // this seems having bad precision
			//}
			result[i] = mat.fullPivLu().solve(rhs);  // find the ith vector

			if (!rhs.isApprox(mat * result[i])) {
				return false;
			}

			if (i + 1 == vector_num)
				break;

			//replace one point with the new vector
			for (SizeType j = 0; j < m_dim; ++j) {
				mat(i + 1, j) = result[i](j, 0);
			}
		}

		return true;
	}

	/* Function to compute the number of sampling such that the error is less
	 * than sigma with probability no less than (1 - 1/e).
	 * @param sigma the bounded error
	 * @param k the dimension to reduce the problem to.
	 */
	mpz_class computeSampleNum(const SizeType sigma, SizeType k) const {
		// dividend = n*...*(n-d+k+1)*(d-k)!
		mpz_class dividend = 1;
		for (SizeType i = m_point_num; i > m_point_num - m_dim + k; --i )
			dividend *= i;
		for (SizeType i = m_dim - k; i > 1; --i)
			dividend *= i;

		// divisor = (sigma+d-k)*...*(sigma+1)*2^{d-k}
		mpz_class divisor;
		mpz_ui_pow_ui(divisor.get_mpz_t(), 2, m_dim - k);
		for (SizeType i = sigma + m_dim - k; i > sigma; --i)
			divisor *= i;

		mpz_class samp_num = dividend / divisor;

		return samp_num;
	}
};

}

#endif /* ABSOLUTEAPX_H_ */
