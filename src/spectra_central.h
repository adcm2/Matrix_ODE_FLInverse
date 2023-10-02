#ifndef SPECTRA_CENTRAL_GUARD_H
#define SPECTRA_CENTRAL_GUARD_H

// #include <math.h>
// #include <omp.h>
// #include <stdio.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <FFTWpp/All>
// #include "./FFTWpp/All"
#include <chrono>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "BlockPreconditioner.h"
#include "blockindex.h"
// #include "filter_header.h"
#include "filter_base.h"
#include "postprocess.h"
#include "postprocessfunctions.h"

namespace modespectrafunctions {

private:
m_rawseismogram, m_filtseismogram;
std::vector<std::complex<double> > m_rawspec, m_filtspec;
bool m_modeinitialized;

Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>
rawspectra(const freq_setup& calcdata, const couplematrix& matdata) {
    // indices
    auto i1 = calcdata.i1();
    auto i2 = calcdata.i2();
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> tmp =
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic,
                      Eigen::Dynamic>::Zero(calcdata.nelem2() + 1,
                                            calcdata.nt() / 2 + 1);
    std::complex<double> myi(0.0, 1.0);
    double oneovertwopi = 1.0 / (2.0 * 3.1415926535);
    std::complex<double> imep =
        static_cast<std::complex<double> >(calcdata.ep());

    // #pragma omp parallel for
    for (int idx = 0; idx < calcdata.nt() / 2 + 1; ++idx) {
        tmp(0, idx) = calcdata.w(idx) * oneovertwopi;   // frequency

        // run through all frequencies and if between f1 and f2 compute
        if (idx > i1 - 1 && idx < i2 + 1) {
            std::complex<double> winp = calcdata.w(idx) - myi * imep;

            // lhs
            Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> vlhs;

            //////////////////////////////////////////////////////////////////////////////////
            // coupling matrix
            Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>
                A(matdata.nelem(), matdata.nelem());

            // declare value of A
            A = matdata.a0() + winp * matdata.a1() + winp * winp * matdata.a2();
            // include diagonal component
            for (int idx = 0; idx < A.rows(); ++idx) {
                A(idx, idx) = A(idx, idx) - winp * winp;
            }

            //////////////////////////////////////////////////////////////////////////////////
            //  rhs and guess
            Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> vrhs, x0;
            vrhs.resize(nelem);
            x0.resize(nelem);

            // find rhs
            for (int idx = 0; idx < nelem; ++idx) {
                vrhs(idx) = matdata.vs(idx) / (myi * winp);
            }

            // finding guess
            Eigen::BlockPreconditioner<std::complex<double> > incond;
            incond.compute(A);
            incond.setblock(A, vecidx[0], vecidx[1]);
            x0 = incond.solve(vrhs);

            //////////////////////////////////////////////////////////////////////////////////

            // declare solver, using diagonal preconditioner for moment
            Eigen::BiCGSTAB<Eigen::Matrix<std::complex<double>, Eigen::Dynamic,
                                          Eigen::Dynamic>,
                            Eigen::BlockPreconditioner<std::complex<double> > >
                solver;

            //////////////////////////////////////////////////////////////////////////////////
            std::vector<int> vecidx;
            vecidx =
                randomfunctions::findindex(winp.real(), wtb, ll, ww.real());
            solver.setTolerance(1.0 * std::pow(10.0, -5));
            solver.compute(A);
            solver.preconditioner().setblock(A, vecidx[0], vecidx[1]);

            //////////////////////////////////////////////////////////////////////////////////

            // solve and return
            vlhs = solver.solveWithGuess(vrhs, x0);

            // find acceleration response using receiver vectors
            tmp.block(1, idx, nelem2, 1) =
                -winp * winp * matdata.vr().transpose() * vlhs;
        };
    };
    return tmp;
};

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
calc_seismogram(const freq_setup& calcdata,
                const Eigen::Matrix<std::complex<double>, Eigen::Dynamic,
                                    Eigen::Dynamic>& vec_rawspec) {
    filterclass::hann freqfilter(calcdata.f1(), calcdata.f2(), 0.1);
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>
        filtspec(vec_rawspec.rows(), vec_rawspec.cols());
    for (int idx = 0; idx < calcdata.nt(); ++idx) {
        for (int idx2 = 0; idx2 < vec_rawspec.rows(); ++idx2) {
            filtspec(idx2, idx) =
                vec_rawspec(idx2, idx) * filters::hannref(calcdata.f(idx),
                                                          calcdata.f1(),
                                                          calcdata.f2(), 0.1);
        }
    }

    // take FT into time domain
};

};   // namespace modespectrafunctions

#endif