#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif
// #ifndef EIGEN_USE_BLAS
// #define EIGEN_USE_BLAS
// #endif
// #ifndef EIGEN_USE_LAPACKE_STRICT
// #define EIGEN_USE_LAPACKE_STRICT
// #endif
#include <math.h>
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
#include "spectra_central.h"

using namespace std::chrono;
int
main() {
    std::string filePath;
    std::string filePath2;
    std::string filePath3;
    std::string pathstring;
    // std::string firstName = "/vector_sr.bin";
    pathstring = std::filesystem::current_path();
    filePath = pathstring + "/matrix.bin";
    filePath2 = pathstring + "/vector_sr.bin";
    filePath3 = pathstring + "/freq_sph.bin";
    auto start = high_resolution_clock::now();
    std::cout << "Hello\n";
    couplematrix mydat(filePath, filePath2, filePath3);
    std::cout << "Hello\n";
    // std::cout << mydat.nelem() << std::endl;
    // std::cout << mydat.a0()(1, 1) << std::endl;
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    std::cout << "Time taken to read in matrices: "
              << duration.count() / 1000000.0 << " seconds" << std::endl;

    start = high_resolution_clock::now();
    // double f1 = 0.1;       // minimum (mHz)
    // double f2 = 1.5;       // maximum (mHz)
    // double dt = 20.0;      // timestep (s)
    // double tout = 512.0;   // time length (hrs)
    // double df0 = 0.05;     // frequency step (mHz)
    // double wtb = 0.05;     // width of target block (mHz)
    // double t1 = 0.0;       // cosine bell start (hrs)
    // double t2 = 512.0;     // cosine bell stop (hrs)
    double f1;
    double f2;
    double dt;
    double tout;
    double df0;
    double wtb;
    double t1;
    double t2;
    double soltol;
    int qex;
    std::cin >> f1 >> f2 >> dt >> tout >> df0 >> wtb >> t1 >> t2 >> soltol >>qex;
    freq_setup myfreq(f1, f2, dt, tout, df0, wtb, t1, t2,qex);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    std::cout << "Time taken to get frequency setup: "
              << duration.count() / 1000000.0 << " seconds" << std::endl;

    // using the new functions for the spectra evaluation
    start = high_resolution_clock::now();
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> rawspec;
    rawspec = modespectrafunctions::rawspectra(myfreq, mydat, soltol);
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    std::cout << "Time taken to do raw calculation: "
              << duration.count() / 1000000.0 << " seconds" << std::endl;

    // finding seismogram
    start = high_resolution_clock::now();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> raw_seis;
    raw_seis =
        modespectrafunctions::calc_seismogram(rawspec, myfreq, mydat.nelem2());
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    std::cout << "Time taken to compute seismogram: "
              << duration.count() / 1000000.0 << " seconds" << std::endl;

    // std::cout << "The value of nt0 is: " << myfreq.nt0() << std::endl;

    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>
        fin_spec;
    start = high_resolution_clock::now();
    fin_spec =
        modespectrafunctions::calc_fspectra(raw_seis, myfreq, mydat.nelem2());
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    std::cout << "Time taken to compute final spectra: "
              << duration.count() / 1000000.0 << " seconds" << std::endl;

    // outputting to files
    std::ofstream myfile;
    std::string outputfilebase = "fspectra.r";
    std::string outputfilename;

    // spectra
    start = high_resolution_clock::now();
    for (int oidx = 0; oidx < mydat.nelem2(); ++oidx) {
        outputfilename = outputfilebase + std::to_string(oidx + 1) + ".out" + ".q" + std::to_string(qex);

        myfile.open(outputfilename, std::ios::trunc);
        for (int idx = myfreq.i12(); idx < myfreq.i22(); ++idx) {
            myfile << std::setprecision(17) << myfreq.f2(idx) * 1000 << ";"
                   << fin_spec(oidx, idx).real() << ";"
                   << fin_spec(oidx, idx).imag() << ";"
                   << std::abs(fin_spec(oidx, idx)) << std::endl;
        }
        myfile.close();
    }

    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    std::cout << "Time taken to output spectra: "
              << duration.count() / 1000000.0 << " seconds" << std::endl;

    // seismogram
    start = high_resolution_clock::now();
    outputfilebase = "tspectra.r";
    for (int oidx = 0; oidx < mydat.nelem2(); ++oidx) {
        outputfilename = outputfilebase + std::to_string(oidx + 1) + ".out";

        myfile.open(outputfilename, std::ios::trunc);
        for (int idx = 0; idx < myfreq.nt(); ++idx) {
            // double tval = idx * mymode.dt / 3600;
            if (myfreq.t(idx) < myfreq.tout()) {
                myfile << std::setprecision(17) << myfreq.t(idx) << ";"
                       << raw_seis(oidx, idx) << std::endl;
            }
        }
        myfile.close();
    }
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    std::cout << "Time taken to output seismograms: "
              << duration.count() / 1000000.0 << " seconds" << std::endl;

    // std::cout << "df0 is: " << myfreq.df0() << std::endl;
    return 0;
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

    // Get duration. Substart timepoints to
    // get duration. To cast it to proper unit
    // use duration cast method

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
    // std::vector<double> myval, myval2;
    // myval.push_back(1.0);
    // myval.push_back(2.0);
    // myval.push_back(2.1);
    // myval2.resize(myval.size());
    // myval2[0] = 1.0;
    // myval2[1] = 2.0;
    // myval2[2] = 1.0;
    // typedef std::vector<double>::iterator ptr;
    // filterclass::hann myhann(0.0, 2.5, 0.1);
    // myhann.filter(myval.begin(), myval.end(), myval2.begin());

    // std::cout << myval2[0] << " " << myval2[1] << " " << myval2[2] <<
    // std::endl;
}
