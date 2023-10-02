// #ifndef EIGEN_DONT_PARALLELIZE
// #define EIGEN_DONT_PARALLELIZE
// #endif
// #ifndef EIGEN_USE_BLAS
// #define EIGEN_USE_BLAS
// #endif
// #ifndef EIGEN_USE_LAPACKE_STRICT
// #define EIGEN_USE_LAPACKE_STRICT
// #endif
#include <omp.h>
#include <stdio.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <vector>

#include "blockindex.h"
#include "filter_base.h"
#include "filter_header.h"
#include "freq_setup.h"
#include "matrix_read.h"
// #include "postprocess.h"
#include "postprocessfunctions.h"

using namespace std::chrono;
int
main() {
    std::vector<double> myval, myval2;
    myval.push_back(1.0);
    myval.push_back(2.0);
    myval.push_back(2.1);
    myval2.resize(myval.size());
    myval2[0] = 1.0;
    myval2[1] = 2.0;
    myval2[2] = 1.0;
    typedef std::vector<double>::iterator ptr;
    filterclass::hann myhann(0.0, 2.5, 0.1);
    myhann.filter(myval.begin(), myval.end(), myval2.begin());

    std::cout << myval2[0] << " " << myval2[1] << " " << myval2[2] << std::endl;

    std::string filePath;
    std::string filePath2;
    std::string filePath3;
    std::string pathstring;
    // std::string firstName = "/vector_sr.bin";
    pathstring = std::filesystem::current_path();
    filePath = pathstring + "/matrix.bin";
    filePath2 = pathstring + "/vector_sr.bin";
    filePath3 = pathstring + "/freq_sph.bin";
    couplematrix mydat(filePath, filePath2, filePath3);
    std::cout << mydat.nelem() << std::endl;
    std::cout << mydat.a0()(1, 1) << std::endl;

    double f1 = 0.1;       // minimum (mHz)
    double f2 = 1.0;       // maximum (mHz)
    double dt = 20.0;      // timestep (s)
    double tout = 512.0;   // time length (hrs)
    double df0 = 0.05;     // frequency step (mHz)
    double wtb = 0.05;     // width of target block (mHz)
    double t1 = 0.0;       // cosine bell start (hrs)
    double t2 = 512.0;     // cosine bell stop (hrs)
    freq_setup myfreq(f1, f2, dt, tout, df0, wtb, t1, t2);

    // int n = 5000;
    // int m;
    // std::cin >> m;
    // Eigen::Matrix<double, Eigen::Dynamic, 1> ll =
    //     Eigen::Matrix<double, Eigen::Dynamic, 1>::Random(n, 1);
    // Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> ww =
    //     Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Random(n, n);
    // Eigen::Matrix<double, Eigen::Dynamic, 1> x;
    // Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> y(n, m);
    // std::vector<int> vecn;
    // for (int i = 0; i < m; ++i) {
    //     vecn.push_back(static_cast<double>(i) / static_cast<double>(m));
    // }

    auto start = high_resolution_clock::now();
    // lusolve.compute(ww);
    // x = lusolve.solve(ll);
    // y = x;
    // printf("Hello from process: %d\n", omp_get_thread_num());
    // Get ending timepoint

    // #pragma omp parallel for
    //     for (int i = 0; i < m; ++i) {
    //         if (i < m) {
    //             Eigen::PartialPivLU<
    //                 Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >
    //                 lusolve;
    //             lusolve.compute(ww * (1 + vecn[i]));
    //             x = lusolve.solve(ll);
    //             y.block(0, i, n, 1) = x;
    //         }
    //         // printf("Hello from process: %d\n", omp_get_thread_num());
    //     }
    auto stop = high_resolution_clock::now();

    // Get duration. Substart timepoints to
    // get duration. To cast it to proper unit
    // use duration cast method
    auto duration = duration_cast<microseconds>(stop - start);

    std::cout << "Time taken by function: " << duration.count() / 1000000.0
              << " seconds" << std::endl;

    return 0;
    // #pragma omp parallel
    //     { printf("Hello from process: %d\n", omp_get_thread_num()); }
    // #pragma omp parallel
    //     { printf("Hello World.\n"); }
    // #pragma omp parallel for ordered schedule(dynamic)
    //     for (int n = 0; n < 10; ++n) {
    // #pragma omp ordered
    //         {
    //             printf(" %d", n);
    //             printf(".\n");
    //         }
    //     }
    // fill out
    // ll.resize(10);
    // ww.resize(10);
    // ll(0) = 0;
    // ll(1) = 1;
    // ll(2) = 1;
    // ll(3) = 1;
    // ll(4) = 2;
    // ll(5) = 2;
    // ll(6) = 2;
    // ll(7) = 2;
    // ll(8) = 2;
    // ll(9) = 2;
    // ww(0) = 0.0;
    // ww(1) = 1.0;
    // ww(2) = 1.0;
    // ww(3) = 2.0;
    // ww(4) = 2.0;
    // ww(5) = 2.0;
    // ww(6) = 3.0;
    // ww(7) = 4.0;
    // ww(8) = 5.0;
    // ww(9) = 6.0;
    // ll.resize(6);
    // ww.resize(6);
    // ll(0) = 2;
    // ll(1) = 2;
    // ll(2) = 1;
    // ll(3) = 3;
    // ll(4) = 3;
    // ll(5) = 4;
    // ww(0) = 0.001943;
    // ww(1) = 0.002382;
    // ww(2) = 0.002538;
    // ww(3) = 0.002944;
    // ww(4) = 0.003683;
    // ww(5) = 0.004066;
    // std::vector<int> validx;
    // validx = randomfunctions::findindex(0.002, 0.0003, ll, ww);
    // std::cout << validx[0] << " " << validx[1] << std::endl;
}
