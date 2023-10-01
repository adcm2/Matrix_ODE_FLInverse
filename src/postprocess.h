#ifndef POSTPROCESS_GUARD_H
#define POSTPROCESS_GUARD_H

#include <FFTWpp>
#include <cassert>
#include <iterator>

// #include "filter_header.h"
#include "filter_base.h"
using namespace filterclass;

namespace postprocess {

enum class ConvertDirection { Forward, Backward };

template <typename INIT>
class PostProcessBase {
    using Float = double;
    using Complex = std::complex<Float>;
    using RealVector = FFTWpp::vector<Float>;
    using ComplexVector = FFTWpp::vector<Complex>;
    using namespace std::complex_literals;

   public:
    /*default constructor*/
    PostProcessBase() : m_isinitialized(false) {}

    /* actual*/
    PostProcessBase(const INIT& vec_inp) : m_isinitialized(true) {
        // length of input and output
        leninp = std::distance(vec_inp.begin(), vec_inp.end());
        // lenout = std::distance(vec_out.begin(), vec_out.end());

        // assign in and out vec
        invec = vec_inp;
        // outvec = vec_out;
    }

    template <typename INIT>
    void transform(ConvertDirection direc) {
        assert(m_isinitialized && "Not initialized");

        // direction of transform
        if (direc == ConvertDirection::Backward) {
            RealVector outFL(2 * (leninp - 1));
            ComplexVector inFL(leninp);
        } else {
            ComplexVector outFL(leninp / 2 + 1);
            RealVector inFL(leninp);
        }

        // Form the plans
        auto flag = FFTWpp::Measure | FFTWpp::Estimate;
        auto inview = FFTWpp::MakeDataView1D(inFL);
        auto outview = FFTWpp::MakeDataView1D(outFL);
        auto fftplan = FFTWpp::Plan(inview, outview, flag);

        // fill out for FT
        auto itinp = invec.begin();
        for (int idx = 0; idx < leninp; ++idx) {
            inFL[idx] = itinp[idx];
        }

        // execute FT
        fftplan.Execute();

        // output
        auto itout = outvec.begin();
        for (int idx = 0; idx < lenout; ++idx) {
            itout[idx] = outFL[idx];
        }

        // return
        return outFL;
    }

   protected:
    mutable bool m_isinitialized;
    auto leninp;
    // auto lenout;
    INIT invec;
    // FIN& outvec;
};

template <typename INIT>
class freq2time : public PostProcessBase<INIT> {
    typedef PostProcessBase<INIT> Base;
    using Base::m_isinitialized;

   public:
    freq2time() : Base() {}

    explicit freq2time(const INIT& vec_inp) : PostProcessBase(vec_inp) {}
}
}   // namespace postprocess
#endif