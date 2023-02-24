#ifndef CEREAL_TYPES_EIGEN_HPP_
#define CEREAL_TYPES_EIGEN_HPP_

#include "cereal/cereal.hpp"
#include <Eigen/Core>

namespace cereal {

    // Serialization for Eigen matrices
    template <class Archive, class _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
    inline void serialize(Archive& archive, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& matrix)
    {
        int rows = matrix.rows();
        int cols = matrix.cols();
        archive(rows, cols);
        matrix.resize(rows, cols);

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                archive(matrix.coeffRef(i, j));
            }
        }
    }

    // Serialization for Eigen vectors
    template <class Archive, class _Scalar, int _Rows, int _Options, int _MaxRows>
    inline void serialize(Archive& archive, Eigen::Matrix<_Scalar, _Rows, 1, _Options, _MaxRows, 1>& vector)
    {
        int rows = vector.rows();
        archive(rows);
        vector.resize(rows);

        for (int i = 0; i < rows; ++i) {
            archive(vector.coeffRef(i));
        }
    }

} // namespace cereal

#endif // CEREAL_TYPES_EIGEN_HPP_
