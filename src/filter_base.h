#ifndef FILTER_BASE_GUARD_H
#define FILTER_BASE_GUARD_H
#include <math.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <chrono>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <cassert>

namespace filterclass {
template <typename xIter>
class filterbase {
    using xtype = std::iter_value_t<xIter>;
   public:
    filterbase() : m_isInitialized(false) {}

    // filterbase(const xtype& p1, const xtype& p2, const xtype& fac)
    //     : m_p1(p1), m_p2(p2), m_fac(fac), m_isInitialized(true) {
    //     m_p11 = p1;
    //     m_p22 = p2;
    //     m_p12 = m_p11 + fac * (m_p22 - m_p11);
    //     m_p21 = m_p22 - fac * (m_p22 - m_p11);
    // };
filterbase(const double& p1, const double& p2, const double& fac)
        : m_p1(p1), m_p2(p2), m_fac(fac), m_isInitialized(true) {
        m_p11 = p1;
        m_p22 = p2;
        m_p12 = m_p11 + fac * (m_p22 - m_p11);
        m_p21 = m_p22 - fac * (m_p22 - m_p11);
    };
    

   protected:
    mutable bool m_isInitialized;
    double m_p1, m_p2, m_fac, m_p11, m_p12, m_p21, m_p22;
    
    bool ltmp;
};

template <typename xIter>
class hann : public filterbase<xIter> {
    typedef filterbase<xIter> Base;
    using filterbase<xIter>::m_fac;
    using filterbase<xIter>::m_isInitialized;
    using filterbase<xIter>::m_p1;
    using filterbase<xIter>::m_p11;
    using filterbase<xIter>::m_p12;
    using filterbase<xIter>::m_p2;
    using filterbase<xIter>::m_p21;
    using filterbase<xIter>::m_p22;
    using Base::ltmp;
    using xtype = std::iter_value_t<xIter>;

   public:
    hann() : Base(){};
    // explicit hann(const xtype &p1, const xtype& p2, const xtype& fac) : Base(p1, p2, fac){};
    explicit hann(const double&p1, const double& p2, const double& fac) : Base(p1, p2, fac){};

    void filter(xIter xbegin, xIter xend, xIter ybegin) {
        assert(m_isInitialized && "Filter not initialized");
        if (xbegin != xend) {
            ltmp = true;
        };
        assert(ltmp && "Not correct inputs");
        xIter it = xbegin;
        int k = 0;
        while (it != xend) {
            ybegin[k] = xbegin[k] * filterindividual(it);
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
    xIter m_xbegin, m_xend, m_ybegin;
    // bool ltmp;
};

}   // namespace filterclass

#endif