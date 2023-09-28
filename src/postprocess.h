#ifndef POSTPROCESS_GUARD_H
#define POSTPROCESS_GUARD_H

#include <FFTWpp>
#include <cassert>

#include "filter_header.h"

namespace postprocess {

template <typename INIT, typename FIN>
class PostProcessBase {
   public:
    /*default constructor*/
    PostProcessBase() : m_isinitialized(false) {}

    template <typename INIT, typename FIN>
    void transform(const INIT& finp, FIN& fout) {
        assert(m_isinitialized && "Not initialized");
        using Float = double;
        using Complex = std::complex<Float>;
        using RealVector = FFTWpp::vector<Float>;
        using ComplexVector = FFTWpp::vector<Complex>;
        using namespace std::complex_literals;

        // set up for FT

        RealVector testFL(nt), checkFL(nt);
        ComplexVector outFL(nt / 2 + 1);

        // Form the plans.
        auto flag = FFTWpp::Measure | FFTWpp::Estimate;

        auto outview = FFTWpp::MakeDataView1D(outFL);
        auto outview2 = FFTWpp::MakeDataView1D(checkFL);
        auto backward_plan = FFTWpp::Plan(outview, outview2, flag);
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> tmpspec;
        Eigen::Matrix<double, 1, Eigen::Dynamic> vecout;

        tmpspec.resize(this->nt / 2 + 1, 1);

        // set up Hann filter
        double fac = 0.1;
        double f11 = this->f1;
        double f22 = this->f2;
        double f12 = f11 + fac * (f22 - f11);
        double f21 = f22 - fac * (f22 - f11);

        // fill out for FT
        for (int idx = 0; idx < nt / 2 + 1; ++idx) {
            // actual frequency
            double finp;
            finp = static_cast<double>(idx) * this->df;
            outFL[idx] =
                rawspec(0, idx) * filters::hann(&finp, &f11, &f12, &f21, &f22);
        }

        // execute FT
        backward_plan.Execute();

        // do corrections
        auto myit2 = checkFL.begin();

        // output
        vecout.resize(1, nt);
        for (int idx = 0; idx < nt; ++idx) {
            double tinp;
            tinp = static_cast<double>(idx) * dt;
            if (tinp < this->tout) {
                vecout(0, idx) = myit2[idx] * exp(this->ep * t[idx]) * df;
            }
        }
        return vecout;
    }

   protected:
    mutable bool m_isinitialized;
};

}   // namespace postprocess
#endif