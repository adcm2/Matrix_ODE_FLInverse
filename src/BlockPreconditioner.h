// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2011-2014 Gael Guennebaud <gael.guennebaud@inria.fr>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef BLOCK_PRECONDITIONERS_H
#define BLOCK_PRECONDITIONERS_H

#include <Eigen/IterativeLinearSolvers>

#include "matrix_replacement.h"
// #include "SolveBlock.h"

namespace Eigen {

template <typename Scalar_>
class BlockPreconditioner {
    typedef Scalar_ Scalar;
    typedef Matrix<Scalar, Dynamic, 1> Vector;

   public:
    typedef typename Vector::StorageIndex StorageIndex;
    enum { ColsAtCompileTime = Dynamic, MaxColsAtCompileTime = Dynamic };

    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    BlockPreconditioner() : m_isInitialized(false) {}

    template <typename MatType>
    explicit BlockPreconditioner(const MatType& mat) : m_invdiag(mat.cols()) {
        compute(mat);
    }

    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////

    EIGEN_CONSTEXPR Index rows() const EIGEN_NOEXCEPT {
        return m_invdiag.size();
    }
    EIGEN_CONSTEXPR Index cols() const EIGEN_NOEXCEPT {
        return m_invdiag.size();
    }

    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////

    // compute call which sets up the preconditioner
    template <typename MatType>
    BlockPreconditioner& compute(const MatType& mat) {
        return factorize(mat);
    }

    // actual call which does something. In particular finds inverse of
    // diagonal entries
    template <typename MatType>
    BlockPreconditioner& factorize(const MatType& mat) {
        // resize inverse of diagonal to correct size
        m_invdiag.resize(mat.cols());

        // go through and find inverse of diagonal
        for (int j = 0; j < mat.outerSize(); ++j) {
            typename MatType::InnerIterator it(mat, j);
            while (it && it.index() != j) ++it;
            if (it && it.index() == j && it.value() != Scalar(0))
                m_invdiag(j) = Scalar(1) / it.value();
            else
                m_invdiag(j) = Scalar(1);
        }

        // set initialized bool to true and return
        m_isInitialized = true;
        return *this;
    }

    // required for the IterativeSolverBase, no meaning here
    template <typename MatType>
    BlockPreconditioner& analyzePattern(const MatType&) {
        return *this;
    }
    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////

    // this is the solve call that is called by the iterative solver base
    // during the iterative solution

    template <typename Rhs, typename Dest>
    void _solve_impl(const Rhs& b, Dest& x) const {
        if (m_isBlockInitialized) {
            for (int j = 0; j < m_invdiag.size(); ++j) {
                if (j < minidx || j > maxidx) {
                    x(j) = b(j) * m_invdiag(j);
                }
            }
            // lu solve
            x.segment(minidx, lenidx) =
                lusolve.solve(b.segment(minidx, lenidx));
        } else {
            // if block for LU decomposition not found/no given
            x = m_invdiag.array() * b.array();
        }
    }
    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    // solve call required by IterativeSolverBase
    template <typename Rhs>
    inline const Solve<BlockPreconditioner, Rhs> solve(
        const MatrixBase<Rhs>& b) const {
        eigen_assert(m_isInitialized &&
                     "BlockPreconditioner is not initialized.");
        eigen_assert(m_invdiag.size() == b.rows() &&
                     "BlockPreconditioner::solve(): invalid number of "
                     "rows of the right hand side matrix b");
        return Solve<BlockPreconditioner, Rhs>(*this, b.derived());
    }

    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    // setting up block LU decomposition
    template <typename MatType>
    void setblock(const MatType& mat, int idx, int idx2) {
        minidx = std::max(std::min(idx, idx2), 0);
        maxidx =
            std::min(std::max(idx, idx2), static_cast<int>(mat.cols()) - 1);
        lenidx = maxidx - minidx + 1;

        if (lenidx > 1) {
            // std::cout << "Hello again\n";
            m_val.resize(lenidx, lenidx);
            // for (int j = 0; j < lenidx; ++j) {
            //     int k = 0;
            //     typename MatType::InnerIterator it(mat, minidx + j);
            //     while (it && it.index() != minidx - 1) ++it;
            //     while (it && it.index() != maxidx) {
            //         ++it;
            //         m_val(j, k) = it.value();
            //         ++k;
            //     }
            // }

            m_val = mat.block(minidx, minidx, lenidx, lenidx);
            // std::cout << m_tmp - m_val << std::endl;
            lusolve.compute(m_val);
            m_isBlockInitialized = true;
        }
    }

    ComputationInfo info() { return Success; }

   protected:
    Vector m_invdiag;
    Matrix<Scalar, Dynamic, Dynamic> m_val;
    int maxidx{0}, minidx{0};
    int lenidx{1};
    bool m_isInitialized;
    bool m_isBlockInitialized{false};
    Eigen::PartialPivLU<Matrix<Scalar, Dynamic, Dynamic> > lusolve;
};

template <typename Scalar_>
class BlockMatFreePreconditioner {
    typedef Scalar_ Scalar;
    typedef Matrix<Scalar, Dynamic, 1> Vector;

   public:
    typedef typename Vector::StorageIndex StorageIndex;
    enum { ColsAtCompileTime = Dynamic, MaxColsAtCompileTime = Dynamic };

    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    BlockMatFreePreconditioner() : m_isInitialized(false) {}

    template <typename MatType>
    explicit BlockMatFreePreconditioner(const MatType& mat)
        : m_invdiag(mat.cols()) {
        compute(mat);
    }

    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////

    EIGEN_CONSTEXPR Index rows() const EIGEN_NOEXCEPT {
        return m_invdiag.size();
    }
    EIGEN_CONSTEXPR Index cols() const EIGEN_NOEXCEPT {
        return m_invdiag.size();
    }

    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////

    // compute call which sets up the preconditioner
    template <typename MatType>
    BlockMatFreePreconditioner& compute(const MatType& mat) {
        return factorize(mat);
    }

    // actual call which does something. In particular finds inverse of
    // diagonal entries
    template <typename MatType>
    BlockMatFreePreconditioner& factorize(const MatType& mat) {
        // resize inverse of diagonal to correct size
        m_invdiag.resize(mat.cols());

        // go through and find inverse of diagonal
        for (int j = 0; j < mat.cols(); ++j) {
            if (mat.coeff(j, j) != Scalar(0))
                m_invdiag(j) = Scalar(1) / mat.coeff(j, j);
            else
                m_invdiag(j) = Scalar(1);
        }

        // set initialized bool to true and return
        m_isInitialized = true;
        return *this;
    }

    // required for the IterativeSolverBase, no meaning here
    template <typename MatType>
    BlockMatFreePreconditioner& analyzePattern(const MatType&) {
        return *this;
    }
    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////

    // this is the solve call that is called by the iterative solver base
    // during the iterative solution

    template <typename Rhs, typename Dest>
    void _solve_impl(const Rhs& b, Dest& x) const {
        if (m_isBlockInitialized) {
            for (int j = 0; j < m_invdiag.size(); ++j) {
                if (j < minidx || j > maxidx) {
                    x(j) = b(j) * m_invdiag(j);
                }
            }
            // lu solve
            x.segment(minidx, lenidx) =
                lusolve.solve(b.segment(minidx, lenidx));
        } else {
            // if block for LU decomposition not found/no given
            x = m_invdiag.array() * b.array();
        }
    }
    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    // solve call required by IterativeSolverBase
    template <typename Rhs>
    inline const Solve<BlockMatFreePreconditioner, Rhs> solve(
        const MatrixBase<Rhs>& b) const {
        eigen_assert(m_isInitialized &&
                     "BlockPreconditioner is not initialized.");
        eigen_assert(m_invdiag.size() == b.rows() &&
                     "BlockPreconditioner::solve(): invalid number of "
                     "rows of the right hand side matrix b");
        return Solve<BlockMatFreePreconditioner, Rhs>(*this, b.derived());
    }

    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    // setting up block LU decomposition
    template <typename MatType>
    void setblock(const MatType& mat, int idx, int idx2) {
        minidx = std::max(std::min(idx, idx2), 0);
        maxidx =
            std::min(std::max(idx, idx2), static_cast<int>(mat.cols()) - 1);
        lenidx = maxidx - minidx + 1;

        if (lenidx > 1) {
            Matrix<Scalar, Dynamic, Dynamic> tmp;
            tmp.resize(lenidx, lenidx);
            for (int i = 0; i < lenidx; ++i) {
                for (int j = 0; j < lenidx; ++j) {
                    tmp(i, j) = mat.coeff(minidx + i, minidx + j);
                }
            }
            lusolve.compute(tmp);
            m_isBlockInitialized = true;
        }
    }

    ComputationInfo info() { return Success; }

   protected:
    Vector m_invdiag;
    // Matrix<Scalar, Dynamic, Dynamic> m_val;
    int maxidx{0}, minidx{0};
    int lenidx{1};
    bool m_isInitialized;
    bool m_isBlockInitialized{false};
    Eigen::PartialPivLU<Matrix<Scalar, Dynamic, Dynamic> > lusolve;
};
}   // namespace Eigen

#endif   // BLOCK_PRECONDITIONERS_H