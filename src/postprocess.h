#ifndef POSTPROCESS_GUARD_H
#define POSTPROCESS_GUARD_H

#include <FFTWpp/All>
#include <cassert>
#include <iterator>

// #include "filter_header.h"
#include "filter_base.h"
using namespace filterclass;

namespace postprocess {
using Float = double;
using Complex = std::complex<Float>;
using RealVector = FFTWpp::vector<Float>;
using ComplexVector = FFTWpp::vector<Complex>;
using namespace std::complex_literals;

enum class ConvertDirection { Forward, Backward };

template <typename INIT, typename FIN>
class PostProcessBase {
   public:
    /*default constructor*/
    PostProcessBase() : m_isinitialized(false) {}

    /* actual*/
    PostProcessBase(const int& len_inp, const int& len_out)
        : m_isinitialized(true), leninp(len_inp), lenout(outlen) {
        inFL.resize(len_inp);
        outFL.resize(len_out);
        // Form the plans
        flag = FFTWpp::Measure | FFTWpp::Estimate;
        inview = FFTWpp::MakeDataView1D(inFL);
        outview = FFTWpp::MakeDataView1D(outFL);
        fftplan = FFTWpp::Plan(inview, outview, flag);
    }

    // template <typename INIT, typename FIN>
    FIN transformcr(const INIT& vecin) {
        assert(m_isinitialized && "Not initialized");
        bool m_rightlen = (vecin.size() == leninp);
        assert(m_rightlen && "Incorrect length");

        // fill out for FT
        auto itinp = vecin.begin();
        for (int idx = 0; idx < leninp; ++idx) {
            inFL[idx] = itinp[idx];
        }

        // execute FT
        fftplan.Execute();

        return this->outFL;
    }

   protected:
    mutable bool m_isinitialized;
    int leninp;
    int lenout;
    FFTWpp::PlanFlag flag;
    FFTWpp::DataView inview;
    FFTWpp::DataView outview;
    FFTWpp::Plan fftplan;
    INIT inFL;
    FIN outFL;
};

// template <typename INIT>
// class freq2time : public PostProcessBase<ComplexVector, RealVector> {
//     typedef PostProcessBase<ComplexVector, RealVector> Base;
//     using Base::invec;
//     // using Base::leninp;
//     using Base::m_isinitialized;

//    public:
//    /*default constructor*/
//     freq2time() : Base() {}

//     freq2time(int len_time) : Base::PostProcessBase(len_time/2+1, len_time) {

//         eiglenout = std::distance(vec_out.begin(), vec_out.end());
//         vec_baseinp.resize(eiglenout);

//         m_lencor = (eiglenout/2 + 1 == leninp);
//         assert(m_lencor && "Incorrect matrix dimensions");
//     }

//     /*defining output*/
//     void transform(const Eigen::Matrix<std::complex<double>, Eigen::Dynamic,
//     Eigen::Dynamic>& vec_in, Eigen::Matrix<double, Eigen::Dynamic,
//     Eigen::Dynamic>& vec_out){
//         auto rownum = vec_in.rows();
//         auto colnum = vec_in.cols();
//         bool m_corrow = (rownum == vec_out.rows());
//         bool m_corcol = (colnum == vec_out.cols()/2+1);
//         assert(m_corrow && "Incorrect rows");
//         assert(m_corcol && "Incorrect columns");

//         RealVector outvec;
//         for (int idx = 0; idx < rownum; ++idx){
//             outvec = Base::transformcr(vec_in.block(idx,0,1,rownum));
//             for (int idx2 = 0; idx2 < colnum; ++idx2){
//                 vec_out(idx,idx2) = outvec[idx2];
//             }

//         }

//         /*time correction*/
//         void correct( Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>&
//         vec_out, const std::vector<double>& vec_time, int tout, double ep,
//         double df){
//             auto colnum = vec_out.cols();
//             for (int idx = 0; idx < rownum; ++idx){

//             }

//         }

//     }

//     private:
//     transform(){
//         RealVector vec_out(2 * )
//     }
//     int eiglenout;
//     bool m_lencor;
//     ComplexVector vec_baseinp;

// }

// template <typename INIT, typename FIN>
// class PostProcessBase {

//    public:
//     /*default constructor*/
//     PostProcessBase() : m_isinitialized(false) {}

//     /* actual*/
//     PostProcessBase(const INIT& vec_inp, const int& outlen) :
//     m_isinitialized(true), lenout(outlen) {
//         // length of input and output
//         leninp = std::distance(vec_inp.begin(), vec_inp.end());

//         // assign in vec
//         invec = vec_inp;
//     }

//     template <typename INIT, typename FIN>
//     void transformcr(FIN& vecfin) {
//         assert(m_isinitialized && "Not initialized");
//         bool m_rightlen = (vecfin.size() == lenout);
//         assert(m_rightlen && "Incorrect length");

//         // direction of transform
//         // INIT inFL(leninp);
//         if (lenout > leninp){
//             RealVector outFL(lenout);
//             ComplexVector inFL(leninp);
//         } else{
//             ComplexVector outFL(lenout);
//             RealVector inFL(leninp);
//         }

//         // Form the plans
//         auto flag = FFTWpp::Measure | FFTWpp::Estimate;
//         auto inview = FFTWpp::MakeDataView1D(inFL);
//         auto outview = FFTWpp::MakeDataView1D(outFL);
//         auto fftplan = FFTWpp::Plan(inview, outview, flag);

//         // fill out for FT
//         auto itinp = invec.begin();
//         for (int idx = 0; idx < leninp; ++idx) {
//             inFL[idx] = itinp[idx];
//         }

//         // execute FT
//         fftplan.Execute();

//         // return
//         auto itout = vecfin.begin();
//         for (int idx = 0; idx < lenout; ++idx){
//             itout[idx] = outFL[idx];
//         }
//         // return outFL;
//     }

//    protected:
//     mutable bool m_isinitialized;
//     auto leninp;
//     auto lenout;
//     INIT invec;
//     // FIN outvec;
// };
}   // namespace postprocess
#endif