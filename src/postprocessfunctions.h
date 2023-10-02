#ifndef POSTPROCESSFUNC_GUARD_H
#define POSTPROCESSFUNC_GUARD_H

#include <FFTWpp>
#include <cassert>
#include <iterator>

// #include "filter_header.h"
#include "filter_base.h"
#include "postprocess.h"
using namespace filterclass;

namespace processfunctions {
using Float = double;
using Complex = std::complex<Float>;
using RealVector = FFTWpp::vector<Float>;
using ComplexVector = FFTWpp::vector<Complex>;
using namespace std::complex_literals;

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
rawfreq2time(const Eigen::Matrix<std::complex<double>, Eigen::Dynamic,
                                 Eigen::Dynamic> &rawspec,
             nt) {   // do Fourier transform

    postprocess::PostProcessBase<ComplexVector, RealVector> myconvert(
        nt / 2 + 1, nt);

    // do corrections
    RealVector vec_out(nt);

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> tmp(rawspec.rows(),
                                                              nt);
    for (int idx = 0; idx < rawspec.rows(); ++idx) {
        vec_out =
            myconvert.transformcr(rawspec.block(idx, 0, 1, rawspec.cols()));
        for (int idx2 = 0; idx2 < nt; ++idx2) {
            tmp(idx, idx2) = vec_out[idx2];
        }
    }

    // return
    return tmp;
};

Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>
rawtime2freq(
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &rawspec,
    nt) {   // do Fourier transform

    postprocess::PostProcessBase<RealVector, ComplexVector> myconvert(
        nt, nt / 2 + 1);

    // do corrections
    ComplexVector vec_out(nt / 2 + 1);

    Eigen::Matrix<std : complex<double>, Eigen::Dynamic, Eigen::Dynamic> tmp(
        rawspec.rows(), nt / 2 + 1);
    for (int idx = 0; idx < rawspec.rows(); ++idx) {
        vec_out =
            myconvert.transformcr(rawspec.block(idx, 0, 1, rawspec.cols()));
        for (int idx2 = 0; idx2 < nt / 2 + 1; ++idx2) {
            tmp(idx, idx2) = vec_out[idx2];
        }
    }

    // return
    return tmp;
};

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
simpfreq2time(const Eigen::Matrix<std::complex<double>, Eigen::Dynamic,
                                  Eigen::Dynamic> &rawspec,
              const std::vector<double> &vec_time, double tout, double ep) {
    int nrow = rawspec.rows();
    int ncol = rawspec.cols();
    int nt = vec_time.size();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> tmp(nrow, nt);
    tmp = rawfreq2time(rawspec, nt);
    for (int idx = 0; idx < nt; ++idx) {
        if (vec_time[idx] < tout) {
            tmp.block(0, idx, nrow, 1) =
                tmp.block(0, idx, nrow, 1) * exp(ep * vec_time[idx]);
        }
    }
};

Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>
simptime2freq(
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &rawspec,
    nt) {   // do Fourier transform

    postprocess::PostProcessBase<RealVector, ComplexVector> myconvert(
        nt, nt / 2 + 1);

    // do corrections
    ComplexVector vec_out(nt);

    Eigen::Matrix<std : complex<double>, Eigen::Dynamic, Eigen::Dynamic> tmp(
        rawspec.rows(), nt / 2 + 1);
    for (int idx = 0; idx < rawspec.rows(); ++idx) {
        vec_out =
            myconvert.transformcr(rawspec.block(idx, 0, 1, rawspec.cols()));
        for (int idx2 = 0; idx2 < nt / 2 + 1; ++idx2) {
            tmp(idx, idx2) = vec_out[idx2];
        }
    }

    // return
    return tmp;
};

}   // namespace processfunctions

#endif