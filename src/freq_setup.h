#ifndef FREQ_SETUP_GUARD_H
#define FREQ_SETUP_GUARD_H

#include <cmath>
#include <iostream>
#include <vector>

class freq_setup {
   public:
    // constructor
    freq_setup(double, double, double, double, double, double, double, double);

    // double functions
    double f(int idx) const { return static_cast<double>(idx) * m_df; }
    double f1() const { return m_f1; }
    double f2() const { return m_f2; }
    double tout() const { return m_tout; }
    double df() const { return m_df; }
    double df0() const { return m_df0; }
    double df2() const { return m_df2; }
    double wtb() const { return m_wtb; }
    double t1() const { return m_t1; }
    double t2() const { return m_t2; }
    double ep() const { return m_ep; }
    double dt() const { return m_dt; }

    // integer functions
    int nt() const { return m_nt; }
    int nt0() const { return m_nt0; }
    int i1() const { return m_i1; }
    int i2() const { return m_i2; }
    int i12() const { return m_i12; }
    int i22() const { return m_i22; }

    // vector returns
    std::vector<double> w() const { return m_w; }
    std::vector<double> t() const { return m_t; }
    double w(int idx) const { return m_w[idx]; }
    double t(int idx) const { return m_t[idx]; }

   private:
    double m_f1, m_f2, m_tout, m_df0, m_wtb, m_t1, m_t2, m_df, m_ep, m_df2,
        m_dt;
    // number of points
    int m_nt, m_nt0, m_i1, m_i2, m_i12, m_i22;

    std::vector<double> m_w;
    std::vector<double> m_t;
};

freq_setup::freq_setup(double f1, double f2, double dt, double tout, double df0,
                       double wtb, double t1, double t2)
    : m_f1{f1 / 1000.0},
      m_f2{f2 / 1000.0},
      m_dt{dt},
      m_df0{df0 / 1000.0},
      m_tout{tout * 3600.0},
      m_t1{t1 * 3600.0},
      m_t2{std::min(t2, tout) * 3600.0},
      m_wtb{wtb * 3.1415926535 / 500.0} {
    // Nyquist frequency
    double fn = 0.5 / dt;

    // check for f2
    if (fn < m_f2) {
        std::cout << "f2 is greater than the Nyquist frequency for the time "
                     "step. Behaviour may be unexpected"
                  << std::endl;
        m_f2 = fn;
    };

    // multipliers for epsilon in inverse FL and q (need to clarify this)
    int mex = 5;
    int qex = 4;
    m_ep = mex / m_tout;
    // David's version
    m_df = m_ep / (6.28318530718 * qex);

    // my version
    // df = 1.0 / this->tout;             // find df
    m_nt = std::ceil(1.0 / (m_df * m_dt));   // find nt
    int ne = static_cast<int>(log(static_cast<double>(m_nt)) / log(2.0) + 1);
    m_nt = pow(2, ne);   // finish increase in density
    std::cout << "nt: " << m_nt << std::endl;
    m_df = 1.0 / (m_nt * m_dt);   // new df

    // finding new values of f1 and f2 and corresponding integer numbers
    m_i1 = std::max(static_cast<int>(std::floor(m_f1 / m_df)), 0);
    m_i2 = static_cast<int>(std::floor(m_f2 / m_df)) + 2;
    m_f1 = (m_i1 - 1) * m_df;
    m_f2 = (m_i2 - 1) * m_df;

    // fill out w, remembering w = 2pi f, and only needing to go up to nt/2+1 as
    // half length of time series
    for (int idx = 0; idx < m_nt / 2 + 1; ++idx) {
        m_w.push_back(2.0 * 3.1415926535 * m_df * static_cast<double>(idx));
    }

    // fill out t
    for (int idx = 0; idx < m_nt; ++idx) {
        m_t.push_back(m_dt * static_cast<double>(idx));
    }
    m_nt0 = floor(1.0 / m_df0 * m_dt);
    if (m_nt0 > m_nt) {
        int ne2 = log(static_cast<double>(m_nt0)) / log(2.0) + 1;
        m_nt0 = pow(2, ne2);
    } else {
        m_nt0 = m_nt;
    }

    m_df2 = 1.0 / (m_nt0 * m_dt);
    m_i12 = std::max(static_cast<int>(floor(m_f1 / m_df2)), 1);
    m_i22 = static_cast<int>(floor(m_f2 / m_df2)) + 1;
};
#endif