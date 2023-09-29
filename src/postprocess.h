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

template <typename INIT, typename FILTVEC, typename FIN>
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
    PostProcessBase(const INIT& vec_inp, const FILTVEC& vec_filt, FIN& vec_out,
                    ConvertDirection direc)
        : m_isinitialized(true) {}
    template <typename INIT, typename FILTVEC, typename FIN>
    void transform(const INIT& vec_inp, const FILTVEC& vec_filt, FIN& vec_out,
                   ConvertDirection direc) {
        assert(m_isinitialized && "Not initialized");

        // length of input and output
        auto leninp = std::distance(vec_inp.begin(), vec_inp.end());
        auto lenout = std::distance(vec_out.begin(), vec_out.end());

        // direction of transform
        if (direc == ConvertDirection::Backward) {
            RealVector outFL(lenout);
            ComplexVector inFL(leninp);
        } else {
            ComplexVector outFL(lenout);
            RealVector inFL(leninp);
        }

        // Form the plans
        auto flag = FFTWpp::Measure | FFTWpp::Estimate;
        auto inview = FFTWpp::MakeDataView1D(inFL);
        auto outview = FFTWpp::MakeDataView1D(outFL);
        auto fftplan = FFTWpp::Plan(inview, outview, flag);

        // fill out for FT
        auto itinp = vec_inp.begin();
        for (int idx = 0; idx < leninp; ++idx) {
            inFL[idx] = itinp[idx];
        }

        // execute FT
        fftplan.Execute();

        // output
        auto itout = vec_out.begin();
        for (int idx = 0; idx < lenout; ++idx) {
            itout[idx] = outFL[idx];
        }
    }

   protected:
    mutable bool m_isinitialized;
    auto leninp;
    auto lenout;
};

}   // namespace postprocess
#endif