#ifndef TYPES_H_
#define TYPES_H_

namespace TD {

typedef std::size_t SizeType;

/*
 * Matrix to store a point set.
 * @tparam DataType type for coordinates of the points.
 */
template<typename DataType>
struct PointSet {
	DataType* m_data;  // pointer to the start of a matrix
	SizeType m_point_num;  // number of points (rows) in the matrix
	SizeType m_dim;  // dimension of the points, i.e. number of columns of the matrix.

	PointSet(DataType* data, SizeType point_num, SizeType dim)
	:m_data(data), m_point_num(point_num),m_dim(dim){}
};

}

#endif /* TYPES_H_ */
