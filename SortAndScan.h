#ifndef SORTANDSCAN_H_
#define SORTANDSCAN_H_

#include <vector>
#include <algorithm>

#include "types.h"

namespace TD {

/*
 * Structure to store a 2D point.
 * @tparam DataType type for coordinates of the points.
 */
template<typename DataType>
struct Point2D {
	DataType x, y;
};

/*
 * SortAndScan is an algorithm for Tukey depth in 2D. Time complexity is bounded
 * by the std::sort. SortAndScan is a template class, so you can use the data type
 * of your choice.
 * @tparam DataType type for coordinates of the points
 */
template<typename DataType>
class SortAndScan {

public:

	/* Function to compute the depth for a point with respect to the other points in the data set.
	 * @param dataset the data set to compute depth with.
	 * @param index the index of the point in the data set (this point is not counted for its depth value).
	 * @return the depth of the point with respect to the other points
	 *  in the data set.
	 */
	SizeType depth(const std::vector<Point2D<DataType>> &dataset, const SizeType index) {
		normalizeDataset(dataset, index);
		return depthOfOrigin();
	}

	/* Operator to compute the depth for a point with respect to the other points in the data set.
	 * @param dataset the data set to compute depth with.
	 * @param index the index of the point in the data set (this point is not counted for its depth value).
	 * @return the depth of the point with respect to the other points
	 *  in the data set.
	 */
	SizeType operator()(const PointSet<DataType> &dataset, const SizeType index) {
		if (dataset.m_dim != 2) {
			throw std::runtime_error("SortAndScan is designated for data in 2d.");
		}
		normalizeDataset(dataset, index);
		return depthOfOrigin();
	}

	/* Function to compute the depth for a point.
	 * @param dataset the data set to compute depth with.
	 * @param p the point to compute depth for.
	 * @return the depth of the point with respect to the data set.
	 */
	SizeType depth(const std::vector<Point2D<DataType>> &dataset,
			const Point2D<DataType> &p) {
		normalizeDataset(dataset, p);
		return depthOfOrigin();
	}

private:
	std::vector<Point2D<DataType>> m_normalized_dataset;

	/* Compare two points pointed to by a and b using the convention that a < b
	 * iff a is hit first when rotating the positive x-axis counterclockwise
	 */
	static bool point2d_lessthan(const Point2D<DataType> &a, const Point2D<DataType> &b) {
		if (b.y == 0 && b.x > 0) {  // b is on the positive x-axis
			return false;
		} else if (a.y == 0 && a.x > 0) {  // only a is on the positive x-axis
			return true;
		} else if (a.y * b.y < 0) {  // a and b are on different sides of the x-aix
			return (a.y > b.y);
		} else {
			return (a.y * b.x - a.x * b.y < 0);  // a < b iff aob is a right turn
		}
	}

	SizeType depthOfOrigin() {
		SizeType num_origins = remove_origin_points();

		SizeType point_num = m_normalized_dataset.size();

		// sort the data
		std::sort(m_normalized_dataset.begin(), m_normalized_dataset.end(), point2d_lessthan);

		// when scanning the data, we will sweep the x-axis counterclockwise
		SizeType num_above = 0;  // number of points above the x-axis

		// Count the number of points that are above and under the x-axis. Points that are on the positive x-axis
		// are counted as above, and the ones on the negative x-axis are counted as below.
		for (SizeType i = 0; i < point_num; i++) {
			if (m_normalized_dataset[i].y > 0)  // above x-axis
				num_above++;
			else if (m_normalized_dataset[i].x > 0 && m_normalized_dataset[i].y == 0)  // on the positive x-axis
				num_above++;
			else
				// counting is done
				break;
		}
		SizeType num_below = point_num - num_above;  // number of points above the x-axis

		// Start to scan
		SizeType count_l = num_above;  // count the number of points on the left of the positive x-axis
		SizeType count_r = num_below;  // count the number of points on the right of the positive x-axis
		SizeType idx_above = 0;
		SizeType idx_below = num_above;
		SizeType depth = std::min(count_l, count_r);  // upper bound of the depth
		while (depth > 0 && (idx_above < num_above || idx_below < point_num)) {
			if (idx_above == num_above) {  // no more points in num_above
				count_r -= point_num - idx_below;  // those in num_below haven't been swept can be ignored for count_r;
				depth = std::min(depth, count_r);  // update depth
				break;
			} else if (idx_below == point_num) {  // no more points in num_below
				count_l -= num_above - idx_above;  // those in num_above haven't been swept can be ignored for count_l;
				depth = std::min(depth, count_l);  // update depth
				break;
			}

			DataType orient = m_normalized_dataset[idx_above].y * m_normalized_dataset[idx_below].x
						- m_normalized_dataset[idx_above].x * m_normalized_dataset[idx_below].y;
			if (orient == 0) {  // colinear
				idx_above++;
				idx_below++;  // don't need to change the counts
			} else if (orient < 0) {  // idx_above,o,idx_below is a right turn
				idx_below++;
				count_r--;
				count_l++;
				depth = std::min(depth, count_r);  // update depth
			} else {  // idx_above,o,idx_below is a left turn
				idx_above++;
				count_r++;
				count_l--;
				depth = std::min(depth, count_l);  // update depth
			}
		}

		m_normalized_dataset.clear();

		return depth + num_origins;  // the duplication is counted for the depth
	}

	void normalizeDataset(const std::vector<Point2D<DataType>> &dataset, const SizeType index) {
		if (index >= dataset.size()) {
			throw std::runtime_error(
				std::string("SortAndScan: Index ") + std::to_string(index) + " is out of bound.");
		}

		// make point at index the origin of the data set
		SizeType point_num = dataset.size() - 1;
		m_normalized_dataset.resize(point_num);
		DataType x = dataset[index].x;
		DataType y = dataset[index].y;

		for (SizeType i = 0; i < index; i++) {
			m_normalized_dataset[i].x = dataset[i].x - x;
			m_normalized_dataset[i].y = dataset[i].y - y;
		}
		for (SizeType i = index + 1; i <= point_num; i++) {
			m_normalized_dataset[i - 1].x = dataset[i].x - x;
			m_normalized_dataset[i - 1].y = dataset[i].y - y;
		}
	}

	void normalizeDataset(const PointSet<DataType> &dataset, const SizeType index) {
		if (index >= dataset.m_point_num) {
			throw std::runtime_error(
				std::string("SortAndScan: Index ") + std::to_string(index) + " is out of bound.");
		}

		const DataType *data = dataset.m_data;

		// make point at index the origin of the data set
		SizeType point_num = dataset.m_point_num - 1;
		m_normalized_dataset.resize(point_num);
		DataType x = *(data + index * 2);
		DataType y = *(data + index * 2 + 1);

		for (SizeType i = 0; i < index; i++) {
			m_normalized_dataset[i].x = *(data + i * 2) - x;
			m_normalized_dataset[i].y = *(data + i * 2 + 1) - y;
		}
		for (SizeType i = index + 1; i <= point_num; i++) {
			m_normalized_dataset[i - 1].x = *(data + i * 2) - x;
			m_normalized_dataset[i - 1].y = *(data + i * 2 + 1) - y;
		}
	}

	void normalizeDataset(const std::vector<Point2D<DataType>> &dataset, const Point2D<DataType> &p) {

		// make p the origin of the data set
		SizeType point_num = dataset.size();
		m_normalized_dataset.resize(point_num);
		DataType x = p.x;
		DataType y = p.y;

		for (SizeType i = 0; i < point_num; i++) {
			m_normalized_dataset[i].x = dataset[i].x - x;
			m_normalized_dataset[i].y = dataset[i].y - y;
		}
	}

	SizeType remove_origin_points() {
		SizeType point_num = m_normalized_dataset.size();
		SizeType i = 0;
		while (i < point_num && (m_normalized_dataset[i].x != 0 || m_normalized_dataset[i].y != 0)) {
			i++;
		}

		for (SizeType j = i; j < point_num; j++) {  // the jth point can be a zero point
			while (j < point_num && m_normalized_dataset[j].x == 0 && m_normalized_dataset[j].y == 0) {  // find the fist non-zero
				j++;
			}
			if (j < point_num) {  // exists
				m_normalized_dataset[i].x = m_normalized_dataset[j].x;
				m_normalized_dataset[i].y = m_normalized_dataset[j].y;
				i++;
			}
		}

		m_normalized_dataset.resize(i);  // remove redundant points.

		return point_num - i;
	}

};

} /* namespace TD */

#endif /* SORTANDSCAN_H_ */
