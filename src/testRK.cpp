#include <math.h>
#include <omp.h>
#include <stdio.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <chrono>
#include <cmath>
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
    // std::cout << "Hello\n";
    couplematrix mydat(filePath, filePath2, filePath3);
    // std::cout << "Hello\n";
    // std::cout << mydat.nelem() << std::endl;
    // std::cout << mydat.a0()(1, 1) << std::endl;
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    // std::cout << "Time taken to read in matrices: "
    //           << duration.count() / 1000000.0 << " seconds" << std::endl;

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
    std::cin >> f1 >> f2 >> dt >> tout >> df0 >> wtb >> t1 >> t2 >> soltol >>
        qex;
    freq_setup myfreq(f1, f2, dt, tout, df0, wtb, t1, t2, qex);

    // how many points
    // std::cout << "Number of points evaluated at is: "
    //           << myfreq.i2() - myfreq.i1() << std::endl;

    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    // std::cout << "Time taken to get frequency setup: "
    //           << duration.count() / 1000000.0 << " seconds" << std::endl;

    // eigensolver
    std::complex<double> myi(0.0, 1.0);
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> fRHS;
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> a0, a1,
        a2, a2inv, Mlarge, vecref, vecf, eiginv, vecr, a2corr;
    start = high_resolution_clock::now();
    a0 = mydat.a0();
    a1 = -myi * mydat.a1();
    a2 = -mydat.a2();
    // std::cout << a0.block(0,0,5,5) << std::endl;
    // std::cout << a1.block(0, 0, 8, 8) << std::endl;
    // std::cout << a2.block(0,0,5,5) << std::endl;
    // #pragma omp parallel
    //     {

    //////////////////////////////////////////////////////////////////////////
    // Form M
    Mlarge.resize(mydat.nelem() * 2, mydat.nelem() * 2);
    std::cout << Mlarge.rows() << " " << Mlarge.cols() << std::endl;

    // top left block
    Mlarge.block(0, 0, mydat.nelem(), mydat.nelem()) =
        Eigen::MatrixXcd::Zero(mydat.nelem(), mydat.nelem());

    // top right block
    Mlarge.block(0, mydat.nelem(), mydat.nelem(), mydat.nelem()) =
        Eigen::MatrixXcd::Zero(mydat.nelem(), mydat.nelem());
    for (int idx = 0; idx < mydat.nelem(); ++idx) {
        // for (int idx2=0; idx2<mydat.nelem(); ++idx2){
        //     Mlarge(idx, idx2 + mydat.nelem()) = 0.0;
        // }
        Mlarge(idx, idx + mydat.nelem()) = 1.0;
    }

    // correct matrix a2
    a2corr = a2;
    for (int idx = 0; idx < mydat.nelem(); ++idx) {
        a2corr(idx, idx) = a2(idx, idx) + 1.0;
    }
    a2inv = a2corr.inverse();

    // bottom left block
    Eigen::MatrixXcd mat_x, mat_v;
    mat_x = -a2inv * a0;
    Mlarge.block(mydat.nelem(), 0, mydat.nelem(), mydat.nelem()) = mat_x;

    // bottom right block
    mat_v = -a2inv * a1;
    Mlarge.block(mydat.nelem(), mydat.nelem(), mydat.nelem(), mydat.nelem()) =
        mat_v;

    // timing
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    std::cout << "Inverse time: " << duration.count() / 1000000.0 << " seconds"
              << std::endl;

    //////////////////////////////////////////////////////////////////////////
    // Form RHS
    fRHS.resize(2 * mydat.nelem());
    fRHS.block(0, 0, mydat.nelem(), 1) = Eigen::VectorXcd::Zero(mydat.nelem());

    Eigen::VectorXcd mat_f;
    mat_f = a2inv * mydat.vs();
    fRHS.block(mydat.nelem(), 0, mydat.nelem(), 1) = mat_f;

    vecr.resize(2 * mydat.nelem(), mydat.vr().cols());

    // top half, ie receiver vector part
    vecr.block(0, 0, mydat.nelem(), mydat.vr().cols()) = mydat.vr();

    // bottom half, ie zeros
    vecr.block(mydat.nelem(), 0, mydat.nelem(), mydat.vr().cols()) =
        Eigen::MatrixXcd::Zero(mydat.nelem(), mydat.vr().cols());
    //////////////////////////////////////////////////////////////////////////
    // eigensolve

    // equiv calc of raw spectra
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> rawspec;
    Eigen::MatrixXd tseis;
    // tseis = mat_tout.real();

    Eigen::MatrixXcd mat_rkout;
    Eigen::VectorXcd vec_k1, vec_k2, vec_k3, vec_k4, vec_xn, vec_xnp;
    vec_xn = Eigen::VectorXcd::Zero(2 * mydat.nelem());
    vec_xnp = Eigen::VectorXcd::Zero(2 * mydat.nelem());
    double val_dt = myfreq.dt();
    int myn = 5000;
    mat_rkout.resize(mydat.nelem2(), myn);
    mat_rkout.block(0, 0, mydat.nelem2(), 1) = vecr * vec_xn;
    std::cout << "Hello" << std::endl;
    start = high_resolution_clock::now();
    for (int idx = 1; idx < myn; ++idx) {
        vec_xn = vec_xnp;
        vec_k1 = Mlarge * vec_xn + fRHS;
        vec_k2 = vec_k1 + 0.5 * dt * Mlarge * vec_k1;
        vec_k3 = vec_k1 + 0.5 * dt * Mlarge * vec_k2;
        vec_k4 = vec_k1 + dt * Mlarge * vec_k3;
        vec_xnp = vec_xn + (1 / 6.0) * dt *
                               (vec_k1 + 2.0 * vec_k2 + 2.0 * vec_k3 + vec_k4);
        // mat_rkout.block(0,idx,mydat.nelem2(),1) = vecr * vec_xnp;
        mat_rkout.block(0, idx, mydat.nelem2(), 1) =
            mydat.vr().transpose() *
            Mlarge.block(mydat.nelem(), 0, mydat.nelem(), 2 * mydat.nelem()) *
            vec_xnp;
        // if (idx%10 == 0) {
        //     // std::cout << "Hello: \n";
        //     std::cout << mat_rkout.block(0,idx,mydat.nelem2(),1) <<
        //     std::endl;
        // }
    }
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    std::cout << "RK time: " << duration.count() / 1000000.0 << " seconds"
              << std::endl;
    // std::cout << tseis.block(0, 0, mydat.nelem2(), 5) << std::endl;
    // std::cout << tseis.block(0, myfreq.nt() - 4, mydat.nelem2(), 4)
    //           << std::endl;

    // std::cout << std::exp(eig_values(0) * static_cast<std::complex<double> >
    // (vec_t[1])) << std::endl; std::cout << "size: " << eig_values.size() <<
    // std::endl;

    // for (int idx = 0; idx < myfreq.nt();++idx){
    //     if (std::isnan(vec_t[idx])){
    //         std::cout << "At idx: " << idx << std::endl;
    //         break;
    //     }
    // }
    // rawspec = modespectrafunctions::rawspectra(myfreq, mydat, soltol);

    // Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>
    // outspec(
    //     mydat.nelem2(), myfreq.nt() / 2 + 1);

    // rawspec = processfunctions::rawtime2freq(tseis, myfreq.nt(),
    // myfreq.dt()); std::cout << rawspec.block(0, 0, mydat.nelem2(), 2) <<
    // std::endl; return outspec;

    // finding seismogram, can use same as for test:
    // Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> raw_seis;
    // raw_seis = processfunctions::filtfreq2time(
    //     rawspec, myfreq.df(), myfreq.f1(), myfreq.f2(), myfreq.nt(), 0.0,
    //     myfreq.dt(), myfreq.tout());
    // raw_seis =
    //     modespectrafunctions::calc_seismogram(rawspec, myfreq,
    //     mydat.nelem2());
    // std::cout << "Hello\n";
    // std::cout << "Hello, it is: " << myfreq.nt0() << std::endl;

    // std::cout << "The matrix P: \n";
    // std::cout << a0.block(0,0,5,5) <<std::endl;
    // time to frequency
    // Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>
    //     fin_spec;
    // fin_spec =
    //     modespectrafunctions::calc_fspectra(raw_seis, myfreq,
    //     mydat.nelem2());
    // // std::cout << "Hello again\n";

    // // outputting to files
    std::ofstream myfile;
    std::string outputfilebase = "fspectra.fullcouple.r";
    std::string outputfilename;

    // spectra
    // for (int oidx = 0; oidx < mydat.nelem2(); ++oidx) {
    //     outputfilename = outputfilebase + std::to_string(oidx + 1) + ".out" +
    //                      ".q" + std::to_string(qex);

    //     myfile.open(outputfilename, std::ios::trunc);
    //     for (int idx = myfreq.i12(); idx < myfreq.i22(); ++idx) {
    //         myfile << std::setprecision(17) << myfreq.f2(idx) * 1000 << ";"
    //                << fin_spec(oidx, idx).real() << ";"
    //                << fin_spec(oidx, idx).imag() << ";"
    //                << std::abs(fin_spec(oidx, idx)) << std::endl;
    //     }
    //     myfile.close();
    // }

    // seismogram
    outputfilebase = "tspectra.RK.r";
    for (int oidx = 0; oidx < mydat.nelem2(); ++oidx) {
        outputfilename = outputfilebase + std::to_string(oidx + 1) + ".out";

        myfile.open(outputfilename, std::ios::trunc);
        for (int idx = 0; idx < myfreq.nt(); ++idx) {
            // double tval = idx * mymode.dt / 3600;
            if (myfreq.t(idx) < myfreq.tout()) {
                myfile << std::setprecision(17) << myfreq.t(idx) << ";"
                       << mat_rkout(oidx, idx) << std::endl;
            }
        }
        myfile.close();
    }

    // }
    return 0;
}