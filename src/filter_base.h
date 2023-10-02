#ifndef FILTER_BASE_GUARD_H
#define FILTER_BASE_GUARD_H
#include <math.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <cassert>
#include <chrono>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

namespace filterclass {
template <typename xtype>
class filterbase {
   public:
    filterbase() : m_isInitialized(false) {}

    filterbase(const xtype& p1, const xtype& p2, const xtype& fac)
        : m_p1(p1), m_p2(p2), m_fac(fac), m_isInitialized(true) {
        m_p11 = p1;
        m_p22 = p2;
        m_p12 = m_p11 + fac * (m_p22 - m_p11);
        m_p21 = m_p22 - fac * (m_p22 - m_p11);
    };

    // filterbase(const double& p1, const double& p2, const double& fac)
    //         : m_p1(p1), m_p2(p2), m_fac(fac), m_isInitialized(true) {
    //         m_p11 = p1;
    //         m_p22 = p2;
    //         m_p12 = m_p11 + fac * (m_p22 - m_p11);
    //         m_p21 = m_p22 - fac * (m_p22 - m_p11);
    //     };

   protected:
    mutable bool m_isInitialized;
    xtype m_p1, m_p2, m_fac, m_p11, m_p12, m_p21, m_p22;

    bool ltmp;
};

template <typename xtype>
class hann : public filterbase<xtype> {
    typedef filterbase<xtype> Base;
    using Base::ltmp;
    using Base::m_fac;
    using Base::m_isInitialized;
    using Base::m_p1;
    using Base::m_p11;
    using Base::m_p12;
    using Base::m_p2;
    using Base::m_p21;
    using Base::m_p22;

   public:
    hann() : Base(){};
    explicit hann(const xtype& p1, const xtype& p2, const xtype& fac)
        : Base(p1, p2, fac){};
    // explicit hann(const double& p1, const double& p2, const double& fac)
    //     : Base(p1, p2, fac){};
    // void filter();

    template <xIter, yIter>
    void filter(xIter xbegin, xIter xend, yIter ybegin) {
        assert(m_isInitialized && "Filter not initialized");
        if (xbegin != xend) {
            ltmp = true;
        };
        assert(ltmp && "Not correct inputs");
        xIter it = xbegin;
        int k = 0;
        while (it != xend) {
            ybegin[k] = ybegin[k] * hann<xtype>::filterindividual(it);
            ++k;
            ++it;
        }
    };

    xtype filterindividual(xIter xpos) {
        if (xpos[0] < m_p11) {
            return static_cast<xtype>(0.0);
        } else if (xpos[0] >= m_p11 && xpos[0] < m_p12) {
            xtype tmp = static_cast<xtype>(3.1415926535) * (xpos[0] - m_p11) /
                        (m_p12 - m_p11);
            return static_cast<xtype>(0.5 * (1.0 - std::cos(tmp)));
        } else if (xpos[0] >= m_p12 && xpos[0] < m_p21) {
            return static_cast<xtype>(1.0);
        } else if (xpos[0] >= m_p21 && xpos[0] < m_p22) {
            xtype tmp = static_cast<xtype>(3.1415926535) * (m_p22 - xpos[0]) /
                        (m_p22 - m_p21);
            return static_cast<xtype>(0.5 * (1.0 - std::cos(tmp)));
        } else {
            return static_cast<xtype>(0.0);
        }
    };

   private:
    // xIter m_xbegin, m_xend;
    // yIter m_ybegin;
    // bool ltmp;
};

}   // namespace filterclass

#endif