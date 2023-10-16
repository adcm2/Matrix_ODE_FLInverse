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

    //how many points
    std::cout << "Number of points evaluated at is: " << myfreq.i2() - myfreq.i1() << std::endl;

    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    std::cout << "Time taken to get frequency setup: "
              << duration.count() / 1000000.0 << " seconds" << std::endl;

    //eigensolver
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> fRHS;
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> a0, a1, a2, a2inv, Mlarge;
    start = high_resolution_clock::now();
    a0 = mydat.a0();
    a1 = mydat.a1();
    a2 = mydat.a2();
// #pragma omp parallel
//     {
    Mlarge.resize(mydat.nelem() * 2,mydat.nelem() * 2 );
    std::cout << Mlarge.rows() << " " << Mlarge.cols() << std::endl;
    Mlarge.block(0,0,mydat.nelem(), mydat.nelem()) = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>::Zero(mydat.nelem(), mydat.nelem());
    for (int idx = 0;idx<mydat.nelem(); ++idx){
        for (int idx2=0; idx2<mydat.nelem(); ++idx2){
            Mlarge(idx, idx2 + mydat.nelem()) = 0.0;
        }
        Mlarge(idx, idx + mydat.nelem()) = 1.0;
    }
    a2inv = a2.inverse();
    Mlarge.block(mydat.nelem(),0,mydat.nelem(), mydat.nelem()) = - a2inv * a0;
    Mlarge.block(mydat.nelem(),mydat.nelem(),mydat.nelem(), mydat.nelem()) = - a2inv * a1;
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    std::cout << "Inverse time: "
              << duration.count() / 1000000.0 << " seconds" << std::endl;
    
    start = high_resolution_clock::now();
    Eigen::ComplexEigenSolver<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> > msolve;
    
    //  Eigen::EigenSolver<Eigen::MatrixXcd> eigensolver;
    msolve.compute(Mlarge);
    // std::cout << Mlarge.rows() << " " << Mlarge.cols() << std::endl;
    Eigen::Vector< std::complex<double>, Eigen::Dynamic>  eig_values = msolve.eigenvalues();
    Eigen::Matrix< std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>  eig_vecs = msolve.eigenvectors();
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    std::cout << "Eigensolution time: "
              << duration.count() / 1000000.0 << " seconds" << std::endl;
    // }
return 0;
}