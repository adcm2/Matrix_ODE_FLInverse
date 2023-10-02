#ifndef FILTER_HEADER_GUARD_H
#define FILTER_HEADER_GUARD_H
#include <math.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <chrono>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace filters {

template <typename xtype>
xtype
hannref(const xtype& f, const xtype& f1, const xtype& f2, const xtype& fac) {
    bool ltmp;
    // double fac = 0.1;
    xtype f11 = f1;
    xtype f22 = f2;
    xtype f12 = f11 + fac * (f22 - f11);
    xtype f21 = f22 - fac * (f22 - f11);
    ltmp = (f11 == 0.0 && f12 == 0.0 && f21 == 0.0 && f22 == 0.0);

    if (!ltmp) {
        if (f < f11) {
            return static_cast<xtype>(0.0);
        } else if (f >= f11 && f < f12) {
            xtype tmp =
                static_cast<xtype>(3.1415926535) * (f - f11) / (f12 - f11);
            return static_cast<xtype>(0.5 * (1.0 - std::cos(tmp)));
        } else if (f >= f12 && f < f21) {
            return static_cast<xtype>(1.0);
        } else if (f >= f21 && f < f22) {
            xtype tmp =
                static_cast<xtype>(3.1415926535) * (f22 - f) / (f22 - f21);
            return static_cast<xtype>(0.5 * (1.0 - std::cos(tmp)));
        } else {
            return static_cast<xtype>(0.0);
        }
    } else {
        return 0.0;
    };
};

// using xvt = std::iter_value_t<xIter>;
template <typename xIter>
std::iter_value_t<xIter>
hann(xIter f, xIter f11, xIter f12, xIter f21, xIter f22) {
    bool ltmp;
    ltmp = (f11[0] == 0.0 && f12[0] == 0.0 && f21[0] == 0.0 && f22[0] == 0.0);
    // std::cout << ltmp << std::endl;
    if (!ltmp) {
        if (f[0] < f11[0]) {
            return static_cast<std::iter_value_t<xIter>>(0.0);
        } else if (f[0] >= f11[0] && f[0] < f12[0]) {
            std::iter_value_t<xIter> tmp =
                static_cast<std::iter_value_t<xIter>>(3.1415926535) *
                (f[0] - f11[0]) / (f12[0] - f11[0]);
            return static_cast<std::iter_value_t<xIter>>(0.5 *
                                                         (1.0 - std::cos(tmp)));
        } else if (f[0] >= f12[0] && f[0] < f21[0]) {
            return static_cast<std::iter_value_t<xIter>>(1.0);
        } else if (f[0] >= f21[0] && f[0] < f22[0]) {
            std::iter_value_t<xIter> tmp =
                static_cast<std::iter_value_t<xIter>>(3.1415926535) *
                (f22[0] - f[0]) / (f22[0] - f21[0]);
            return static_cast<std::iter_value_t<xIter>>(0.5 *
                                                         (1.0 - std::cos(tmp)));
        } else {
            return static_cast<std::iter_value_t<xIter>>(0.0);
        }
    } else {
        return 0.0;
    };
};

}   // namespace filters

#endif