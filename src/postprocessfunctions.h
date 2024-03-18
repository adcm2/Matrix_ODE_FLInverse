#ifndef POSTPROCESSFUNC_GUARD_H
#define POSTPROCESSFUNC_GUARD_H

#include <FFTWpp/Ranges>
#include <cassert>
#include <iterator>
#include <ranges>

// #include "filter_header.h"
#include "filter_base.h"
// #include "postprocess.h"
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
             const int nt) {   // do Fourier transform

    // do corrections
    using namespace FFTWpp;

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> tmp(rawspec.rows(),
                                                              nt);
    ComplexVector inFL(nt / 2 + 1);
    RealVector outFL(nt);

    // Form the plans
    // auto flag = FFTWpp::Measure | FFTWpp::Estimate;
    // auto inview = FFTWpp::MakeDataView1D(inFL);
    // auto outview = FFTWpp::MakeDataView1D(outFL);
    // auto fftplan = FFTWpp::Plan(inview, outview, flag);
    auto planForward = Ranges::Plan(Ranges::View(inFL), Ranges::View(outFL), FFTWpp::Measure);

    for (int idx = 0; idx < rawspec.rows(); ++idx) {
        // auto itinp = .begin();
        for (int idx2 = 0; idx2 < nt / 2 + 1; ++idx2) {
            inFL[idx2] = rawspec(idx, idx2);
        }

        // do FFT
        // fftplan.Execute();
        planForward.Execute();

        for (int idx2 = 0; idx2 < nt; ++idx2) {
            tmp(idx, idx2) = outFL[idx2];
        }
    }

    // return
    return tmp;
};

Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>
rawtime2freq(
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &rawspec,
    const int &nt, double dt) {   // do Fourier transform

using namespace FFTWpp;
    // for FFT
    RealVector inFL(nt);
    ComplexVector outFL(nt / 2 + 1);

    // Form the plans
    // auto flag = FFTWpp::Measure | FFTWpp::Estimate;
    // auto inview = FFTWpp::MakeDataView1D(inFL);
    // auto outview = FFTWpp::MakeDataView1D(outFL);
    // auto fftplan = FFTWpp::Plan(inview, outview, flag);
    auto planForward = Ranges::Plan(Ranges::View(inFL), Ranges::View(outFL), FFTWpp::Measure);

    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> tmp(
        rawspec.rows(), nt / 2 + 1);
    for (int idx = 0; idx < rawspec.rows(); ++idx) {
        for (int idx2 = 0; idx2 < nt; ++idx2) {
            inFL[idx2] = rawspec(idx, idx2);
        }
        // fftplan.Execute();
        planForward.Execute();
        for (int idx2 = 0; idx2 < nt / 2 + 1; ++idx2) {
            tmp(idx, idx2) = outFL[idx2] * dt;
        }
    }

    // return
    return tmp;
};

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
filtfreq2time(const Eigen::Matrix<std::complex<double>, Eigen::Dynamic,
                                  Eigen::Dynamic> &rawspec,
              const double df, const double f1, const double f2, const int nt,
              const double ep, const double dt, const double tout) {
    // size of matrices
    int nrow = rawspec.rows();
    int ncol = rawspec.cols();

    // declaring temporaries
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> tmp(nrow, nt);
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> tmpraw;
    tmpraw = rawspec;

    // filter raw spectrum
    for (int idx = 0; idx < nt / 2 + 1; ++idx) {
        tmpraw.block(0, idx, nrow, 1) *=
            filters::hannref(df * idx, f1, f2, 0.1);
    }

    // do FFT
    tmp = rawfreq2time(tmpraw, nt);

    // undo effect of frequency shift
    for (int idx = 0; idx < nt; ++idx) {
        if (dt * idx < tout) {
            tmp.block(0, idx, nrow, 1) *= exp(ep * dt * idx) * df;
        }
    }
    return tmp;
};

Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>
simptime2freq(
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &rawspec,
    const double dt, const double t2) {   // do Fourier transform

    // declarations
    int nrow = rawspec.rows();
    int nt = rawspec.cols();
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> tmp(
        rawspec.rows(), nt / 2 + 1);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> tmpraw;

    // filter
    tmpraw = rawspec;
    for (int idx = 0; idx < nt; ++idx) {
        tmpraw.block(0, idx, nrow, 1) *=
            filters::hannref(idx * dt, 0.0, t2, 0.5);
    }

    // do conversion
    tmp = rawtime2freq(tmpraw, nt, dt);

    // return
    return tmp;
};

Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>
fulltime2freq(
    const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> &rawspec,
    const freq_setup &calcdata) {   // do Fourier transform

    // declarations
    int nrow = rawspec.rows();
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> tmp;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> tmpraw =
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(
            nrow, calcdata.nt0());

    // filter
    tmpraw.block(0, 0, nrow, calcdata.nt()) = rawspec;

    // do conversion
    tmp = simptime2freq(tmpraw, calcdata.dt(), calcdata.t2());

    // return
    return tmp;
};

}   // namespace processfunctions

#endif