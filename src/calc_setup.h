#ifndef FREQ_HEADER_GUARD_H
#define FREQ_HEADER_GUARD_H

#include <math.h>

#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <chrono>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "matrix_header.h"

class freq_setup {
   public:
    // constructor
    freq_setup(double, double, double, double, double, double, double, double);

    // operator
    double operator()(double);

    // frequency of particular point
    // std::complex<double> freq_value(int);
    // Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>
    // fspectra();

    // number of points
    int nt;
    int i1, i2;
    double df, ep;

   private:
    double f1, f2, dt, tout, df0, wtb, t1, t2;
};

freq_setup::freq_setup(double f1, double f2, double dt, double tout, double df0,
                       double wtb, double t1, double t2)
    : f1{f1 / 1000.0},
      f2{f2 / 1000.0},
      dt{dt},
      df0{df0 / 1000.0},
      tout{tout * 3600.0},
      t1{t1 * 3600.0},
      t2{std::min(t2, tout) * 3600.0},
      wtb{wtb * 3.1415926535 / 500.0} {
    double fn = 0.5 / dt;
    if (fn < f2) {
        std::cout << "f2 is greater than the Nyquist frequency for the time "
                     "step. Behaviour may be unexpected"
                  << std::endl;
        this->f2 = fn;
    };
    int mex = 5;
    int qex = 4;

    ep = mex / freq_setup::tout;

    df = ep / (6.28318530718 * qex);
    // std::cout << "df: " << df << std::endl;
    // std::cout << "ep: " << ep << std::endl;
    // std::cout << "dt: " << dt << std::endl;
    // std::cout << "t/dt: " << freq_setup::tout / dt << std::endl;

    nt = std::ceil(1.0 / (df * dt));
    std::cout << "nt: " << nt << std::endl;
    int ne = static_cast<int>(log(static_cast<double>(nt)) / log(2.0) + 1);
    std::cout << "ne: " << ne << std::endl;
    nt = pow(2, ne);

    df = 1.0 / (nt * dt);

    i1 = std::max(static_cast<int>(std::floor(freq_setup::f1 / df)), 2);
    freq_setup::f1 = (i1 - 1) * df;
    std::cout << "nt: " << nt << std::endl;
    std::cout << "f1: " << f1 << std::endl;
    i2 = static_cast<int>(std::floor(freq_setup::f2 / df)) + 2;

    freq_setup::f2 = (i2 - 1) * df;
    std::cout << freq_setup::f1 << " " << freq_setup::f2 << std::endl;
};

double
freq_setup::operator()(double x) {
    //  std::cout << f1 << std::endl;
    //  std::cout << f2 << std::endl;
    //  std::cout << dt << std::endl;
    //  std::cout << tout << std::endl;
    //  std::cout << df0 << std::endl;
    //  std::cout << wtb << std::endl;
    //  std::cout << t1 << std::endl;
    //  std::cout << t2 << std::endl;
    return x;
};

// Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>
//         freq_setup::fspectra(){

//         };

// std::complex<double>
// freq_setup::freq_value(int n) {
//     using namespace std::complex_literals;
//     // std::cout << i1 << " " << i2 << std::endl;
//     if (n < i1) {
//         return
//     };
//     if (n > i2) {
//         std::cout << "n is too large\n" << std::endl;
//     };
//     return n * df + ep * 1i;
// };
#endif