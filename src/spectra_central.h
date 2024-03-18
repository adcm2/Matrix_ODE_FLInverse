#ifndef SPECTRA_CENTRAL_GUARD_H
#define SPECTRA_CENTRAL_GUARD_H

#include <math.h>
#include <omp.h>
#include <stdio.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <FFTWpp/Ranges>
#include <complex>
#include <vector>

#include "BlockPreconditioner.h"
#include "blockindex.h"
#include "postprocessfunctions.h"

namespace modespectrafunctions {
using complexd = std::complex<double>;

Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>
rawspectra(const freq_setup& calcdata, const couplematrix& matdata,
           const double soltol) {
    // indices
    auto i1 = calcdata.i1();
    auto i2 = calcdata.i2();
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> tmp =
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic,
                      Eigen::Dynamic>::Zero(matdata.nelem2(),
                                            calcdata.nt() / 2 + 1);

    //////////////////////////////////////////////////////////////////////////////////
    // some simple values
    std::complex<double> myi(0.0, 1.0);
    double oneovertwopi = 1.0 / (2.0 * 3.1415926535897932);
    std::complex<double> imep =
        static_cast<std::complex<double> >(calcdata.ep());
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> A(
        matdata.nelem(), matdata.nelem());
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> mat_a0(
        matdata.nelem(), matdata.nelem()),
        mat_a1(matdata.nelem(), matdata.nelem()),
        mat_a2(matdata.nelem(), matdata.nelem());
    mat_a0 = matdata.a0();
    mat_a1 = matdata.a1();
    mat_a2 = matdata.a2();
    // shared(mat_a0, mat_a1, mat_a2)
//////////////////////////////////////////////////////////////////////////////////
#pragma omp parallel private(A) shared(mat_a0, mat_a1, mat_a2)
    {
#pragma omp for schedule(dynamic, 10)
        for (int idx = i1; idx < i2; ++idx) {
            // coupling matrix

            // tmp(0, idx) = calcdata.w(idx) * oneovertwopi;   // frequency

            // run through all frequencies and if between f1 and f2 compute
            // if (idx > i1  && idx < i2) {
            std::complex<double> winp = calcdata.w(idx) - myi * imep;

            // lhs
            Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> vlhs;

            //////////////////////////////////////////////////////////////////////////////////

            // declare value of A
            // A = matdata.a0() + winp * matdata.a1() + winp * winp *
            // matdata.a2();
            A = mat_a0 + winp * mat_a1 + winp * winp * mat_a2;

            // include diagonal component
            for (int idx = 0; idx < A.rows(); ++idx) {
                A(idx, idx) = A(idx, idx) - winp * winp;
            }

            //////////////////////////////////////////////////////////////////////////////////
            //  rhs and guess
            Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> vrhs, x0;
            vrhs.resize(matdata.nelem());
            x0.resize(matdata.nelem());

            // find rhs
            for (int idx = 0; idx < matdata.nelem(); ++idx) {
                vrhs(idx) = matdata.vs(idx) / (myi * winp);
            }

            // finding guess
            std::vector<int> vecidx;
            vecidx = randomfunctions::findindex(
                winp.real(), calcdata.wtb(), matdata.ll(), matdata.ww().real());
            Eigen::BlockPreconditioner<std::complex<double> > incond;
            incond.compute(A);
            incond.setblock(A, vecidx[0], vecidx[1]);
            x0 = incond.solve(vrhs);

            // LU decomposition
            //  Eigen::FullPivLU<Eigen::Matrix<std::complex<double>,
            //  Eigen::Dynamic, Eigen::Dynamic> > fsolve; fsolve.compute(A);
            //  vlhs = fsolve.solve(vrhs);

            //////////////////////////////////////////////////////////////////////////////////

            // declare solver, using diagonal preconditioner for moment
            Eigen::BiCGSTAB<Eigen::Matrix<std::complex<double>, Eigen::Dynamic,
                                          Eigen::Dynamic>,
                            Eigen::BlockPreconditioner<std::complex<double> > >
                solver;

            //////////////////////////////////////////////////////////////////////////////////

            solver.setTolerance(soltol);
            solver.compute(A);
            solver.preconditioner().setblock(A, vecidx[0], vecidx[1]);

            //////////////////////////////////////////////////////////////////////////////////

            // solve and return
            vlhs = solver.solveWithGuess(vrhs, x0);

            // find acceleration response using receiver vectors
            tmp.block(0, idx, matdata.nelem2(), 1) =
                -winp * winp * matdata.vr().transpose() * vlhs;
            // };
        };
    };
    return tmp;
};

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
calc_seismogram(const Eigen::Matrix<std::complex<double>, Eigen::Dynamic,
                                    Eigen::Dynamic>& vec_rawspec,
                const freq_setup& calcdata, const int nelem2) {
    // declare filtered matrix
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>
        filtspec(nelem2, calcdata.nt() / 2 + 1);

    // take FT into time domain
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> seisout(
        nelem2, calcdata.nt());

    seisout = processfunctions::filtfreq2time(
        vec_rawspec, calcdata.df(), calcdata.f1(), calcdata.f2(), calcdata.nt(),
        calcdata.ep(), calcdata.dt(), calcdata.tout());

    return seisout;
};

Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>
calc_fspectra(
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& seisout,
    const freq_setup& calcdata, const int nelem2) {
    // declaration
    double twopi = 2.0 * 3.1415926535897932;
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> outspec(
        nelem2, calcdata.nt0() / 2 + 1);

    outspec = processfunctions::fulltime2freq(seisout, calcdata);
    return outspec;
}

};   // namespace modespectrafunctions

#endif