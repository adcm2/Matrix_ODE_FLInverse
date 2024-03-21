#ifndef MATRIX_MULT_REPLACE_GUARD_H
#define MATRIX_MULT_REPLACE_GUARD_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <cassert>
#include <iostream>
#include <unsupported/Eigen/IterativeSolvers>

template <typename MRScalar>
class MatrixReplaceFT;
using Eigen::SparseMatrix;

namespace Eigen {
namespace internal {
// MatrixReplaceFT looks-like a SparseMatrix, so let's inherits its traits:
template <typename MRScalar>
struct traits<MatrixReplaceFT<MRScalar>>
    : public Eigen::internal::traits<Eigen::SparseMatrix<MRScalar>> {};
}   // namespace internal
}   // namespace Eigen

// Example of a matrix-free wrapper from a user type to Eigen's compatible type
// For the sake of simplicity, this example simply wrap a Eigen::SparseMatrix.
template <typename MRScalar>
class MatrixReplaceFT : public Eigen::EigenBase<MatrixReplaceFT<MRScalar>> {
   public:
    // Required typedefs, constants, and method:
    using Scalar = MRScalar;
    using RealScalar = double;
    using StorageIndex = int;
    using Index = Eigen::EigenBase<MatrixReplaceFT<MRScalar>>::Index;
    // using namespace GSHTrans;
    // using Real = double;

    using MATRIX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

    enum {
        ColsAtCompileTime = Eigen::Dynamic,
        MaxColsAtCompileTime = Eigen::Dynamic,
        IsRowMajor = false
    };

    Index rows() const { return _mat_a0->rows(); }
    Index cols() const { return _mat_a0->cols(); }
    // Index rows() const { return this->_matlen * this->_size0; }
    // Index cols() const { return this->_matlen * this->_size0; }

    template <typename Rhs>
    Eigen::Product<MatrixReplaceFT<Scalar>, Rhs, Eigen::AliasFreeProduct>
    operator*(const Eigen::MatrixBase<Rhs> &x) const {
        return Eigen::Product<MatrixReplaceFT<Scalar>, Rhs,
                              Eigen::AliasFreeProduct>(*this, x.derived());
    }

    // Custom API:
    MatrixReplaceFT() : _mat_a0(0), _mat_a1(0), _mat_a2(0) {}

    void attachMatrices(const MATRIX &a0, const MATRIX &a1, const MATRIX &a2) {
        assert((a0.rows() == a1.rows() && a0.rows() == a2.rows()) &&
               "Different matrix sizes");
        _mat_a0 = &a0;
        _mat_a1 = &a1;
        _mat_a2 = &a2;
    }
    // return actual matrices
    const MATRIX a0() const { return *_mat_a0; }
    const MATRIX a1() const { return *_mat_a1; }
    const MATRIX a2() const { return *_mat_a2; }

    // return references to matrices
    const MATRIX *a0p() const { return _mat_a0; }
    const MATRIX *a1p() const { return _mat_a1; }
    const MATRIX *a2p() const { return _mat_a2; }

    // return value of matrix at particular index
    const MRScalar coeff(int idxi, int idxj) const {
        assert((idxi < a0.rows() && idxj < a0.cols()) && "Incorret index");
        MRScalar tmp;

        if (idxi == idxj) {
            tmp -= this->w() * this->w();
        }
        tmp += (this->a0p())->coeff(idxi, idxj) +
               (this->a1p())->coeff(idxi, idxj) * this->w() +
               (this->a2p())->coeff(idxi, idxj) * this->w() * this->w();
        return tmp;
    }
    void frequency(MRScalar winp) { _w = winp; };
    MRScalar w() const { return _w; }

   private:
    const SparseMatrix<Scalar> *mp_mat;
    const MATRIX *_mat_a0, *_mat_a1, *_mat_a2;

    MRScalar _w;
};
namespace Eigen {
namespace internal {

template <typename MRScalar, typename Rhs>
struct generic_product_impl<MatrixReplaceFT<MRScalar>, Rhs, SparseShape,
                            DenseShape,
                            GemvProduct>   // GEMV stands for matrix-vector
    : generic_product_impl_base<
          MatrixReplaceFT<MRScalar>, Rhs,
          generic_product_impl<MatrixReplaceFT<MRScalar>, Rhs>> {
    // typedef typename Product<MatrixReplacement<MRScalar>, Rhs>::Scalar
    // Scalar;
    using Scalar = Product<MatrixReplaceFT<MRScalar>, Rhs>::Scalar;
    using Real = Product<MatrixReplaceFT<MRScalar>, Rhs>::RealScalar;
    using Complex = std::complex<Real>;
    // using MATRIX = MatrixReplaceFT::MATRIX;

    // if we have a vector x, we are defining x + alpha * lhs * rhs
    template <typename Dest>
    static void scaleAndAddTo(Dest &dst, const MatrixReplaceFT<MRScalar> &lhs,
                              const Rhs &rhs, const Scalar &alpha) {
        // This method should implement "dst += alpha * lhs * rhs" inplace,
        // however, for iterative solvers, alpha is always equal to 1, so let's
        // not bother about it.
        assert(alpha == Scalar(1) && "scaling is not implemented");
        EIGEN_ONLY_USED_FOR_DEBUG(alpha);
        // Eigen::MatrixXcd A =
        //     lhs.a0() + lhs.w() * lhs.a1() + std::pow(lhs.w(), 2.0) *
        //     lhs.a2();
        // Here we could simply call dst.noalias() += lhs.my_matrix() * rhs,
        // but let's do something fancier (and less efficient):
        // dst.noalias() += lhs.a0() * rhs;

        // dst.noalias() += alpha *
        //                      (*(lhs.a0p()) + lhs.w() * (*(lhs.a1p())) +
        //                       lhs.w() * lhs.w() * (*(lhs.a2p()))) *
        //                      rhs -
        //                  alpha * lhs.w() * lhs.w() * rhs;

        // dst.noalias() += alpha * (*(lhs.a0p())) * rhs +
        //                  alpha * lhs.w() * (*(lhs.a1p())) * rhs +
        //                  alpha * lhs.w() * lhs.w() * (*(lhs.a2p())) * rhs -
        //                  alpha * lhs.w() * lhs.w() * rhs;

        dst.noalias() +=
            alpha * ((*(lhs.a0p())) * rhs + lhs.w() * (*(lhs.a1p())) * rhs +
                     lhs.w() * lhs.w() * (*(lhs.a2p())) * rhs -
                     lhs.w() * lhs.w() * rhs);
        // dst.noalias() += alpha *
        //                  ((lhs.a0()) + lhs.w() * ((lhs.a1())) +
        //                   lhs.w() * lhs.w() * ((lhs.a2()))) *
        //                  rhs;
        // dst.noalias() += alpha * A * rhs;
        // dst.noalias() -= alpha * lhs.w() * lhs.w() * rhs;
    }
};
}   // namespace internal
}   // namespace Eigen

#endif