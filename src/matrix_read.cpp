#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

// #ifndef EIGEN_USE_BLAS
// #define EIGEN_USE_BLAS
// #endif
#include <math.h>
#include <omp.h>
#include <stdio.h>

#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <chrono>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// #include "calc_setup.h"
#include "matrix_header.h"
// #include "spectra_header.h"
#include "spectra_eval.h"
#include "matrix_read.h"

using std::cout;
using std::endl;
using std::ifstream;
using std::ios;
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

    // std::cout << pathtest2 << std::endl;

    // modematrix testclass(filePath);
    // double f1 = 0.1;       // minimum (mHz)
    // double f2 = 1.0;       // maximum (mHz)
    // double dt = 20.0;      // timestep (s)
    // double tout = 512.0;   // time length (hrs)
    // double df0 = 0.05;     // frequency step (mHz)
    // double wtb = 0.05;      // width of target block (mHz)
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
    std::cin >> f1 >> f2 >> dt >> tout >> df0 >> wtb >> t1 >> t2;
    // freq_setup mytest(f1, f2, dt, tout, df0, wtb, t1, t2);
    // mytest(5);
    // Get starting timepoint
    
    auto start = high_resolution_clock::now();
    modespectra mymode(filePath, filePath2, filePath3, f1, f2, dt, tout, df0,
                       wtb, t1, t2);
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> mytmp;
    mytmp = mymode.fspectra();

    // outputting to files
    std::ofstream myfile;
    std::string outputfilebase = "fspectra.r";
    std::string outputfilename;
    for (int oidx = 0; oidx < mymode.nelem2; ++oidx) {
        outputfilename = outputfilebase + std::to_string(oidx + 1) + ".out";

        myfile.open(outputfilename, std::ios::trunc);
        for (int idx = mymode.i12; idx < mymode.i22; ++idx) {
            myfile << mymode.df2 * idx * 1000.0 << ";"
                   << mymode.fspec(oidx, idx).real() << ";"
                   << mymode.fspec(oidx, idx).imag() << ";"
                   << std::abs(mymode.fspec(oidx, idx)) << std::endl;
        }
        myfile.close();
    }
    outputfilebase = "tspectra.r";
    for (int oidx = 0; oidx < mymode.nelem2; ++oidx) {
        outputfilename = outputfilebase + std::to_string(oidx + 1) + ".out";

        myfile.open(outputfilename, std::ios::trunc);
        for (int idx = 0; idx < mymode.nt; ++idx) {
            double tval = idx * mymode.dt / 3600;
            // if (tval < tout){
            myfile << mymode.dt * idx << ";" << mymode.tseis(oidx, idx)
                   << std::endl;
        }
        myfile.close();
    }

    // Get ending timepoint
    auto stop = high_resolution_clock::now();

    // Get duration. Substart timepoints to
    // get duration. To cast it to proper unit
    // use duration cast method
    auto duration = duration_cast<microseconds>(stop - start);

    std::cout << "Time taken by function: " << duration.count() / 1000000.0
              << " seconds" << std::endl;
    return 0;
    // test spectra
    // testclass.fspectra(mytest.freq_value(1));
}