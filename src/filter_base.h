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
template <typename xIter, typename yIter>
class filterbase {
    using xtype = std::iter_value_t<xIter>;
    using ytype = std::iter_value_t<yIter>;

   public:
    filterbase() : m_isInitialized(false) {}

    filterbase(xIter xbegin, xIter xend, yIter ybegin, const xtype& p1,
               const xtype& p2, const xtype& fac)
        : m_xbegin(xbegin),
          m_xend(xend),
          m_ybegin(ybegin),
          m_p1(p1),
          m_p2(p2),
          m_fac(fac),
          m_isInitialized(true) {
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
    xIter m_xbegin, m_xend;
    yIter m_ybegin;
    bool ltmp;
};

template <typename xIter, typename yIter>
class hann : public filterbase<xIter, yIter> {
    typedef filterbase<xIter, yIter> Base;
    using Base::ltmp;
    using Base::m_fac;
    using Base::m_isInitialized;
    using Base::m_p1;
    using Base::m_p11;
    using Base::m_p12;
    using Base::m_p2;
    using Base::m_p21;
    using Base::m_p22;
    using Base::m_xbegin;
    using Base::m_xend;
    using Base::m_ybegin;
    using xtype = std::iter_value_t<xIter>;

   public:
    hann() : Base(){};
    explicit hann(xIter xbegin, xIter xend, yIter ybegin, const xtype& p1,
                  const xtype& p2, const xtype& fac)
        : Base(xbegin, xend, ybegin, p1, p2, fac){};
    // explicit hann(const double& p1, const double& p2, const double& fac)
    //     : Base(p1, p2, fac){};
    // void filter();

    void filter() {
        assert(m_isInitialized && "Filter not initialized");
        if (m_xbegin != m_xend) {
            ltmp = true;
        };
        assert(ltmp && "Not correct inputs");
        xIter it = m_xbegin;
        int k = 0;
        while (it != m_xend) {
            m_ybegin[k] =
                m_ybegin[k] * hann<xIter, yIter>::filterindividual(it);
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