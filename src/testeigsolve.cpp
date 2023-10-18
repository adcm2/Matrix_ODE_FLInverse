#include <math.h>
#include <omp.h>
#include <stdio.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <vector>
#include <cmath>

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
    std::complex<double> myi(0.0, 1.0);
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> fRHS;
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> a0, a1, a2, a2inv, Mlarge, vecref, vecf, eiginv, vecr, a2corr;
    start = high_resolution_clock::now();
    a0 = mydat.a0();
    a1 = myi * mydat.a1();
    a2 = -mydat.a2();
    // std::cout << a0.block(0,0,5,5) << std::endl;
    std::cout << a1.block(0,0,8,8) << std::endl;
    // std::cout << a2.block(0,0,5,5) << std::endl;
// #pragma omp parallel
//     {

    ////////////////////////////////////////////////////////////////////////// 
    // Form M
    Mlarge.resize(mydat.nelem() * 2,mydat.nelem() * 2 );
    std::cout << Mlarge.rows() << " " << Mlarge.cols() << std::endl;
    Mlarge.block(0,0,mydat.nelem(), mydat.nelem()) = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>::Zero(mydat.nelem(), mydat.nelem());
    Mlarge.block(0,mydat.nelem(),mydat.nelem(), mydat.nelem()) = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>::Zero(mydat.nelem(), mydat.nelem());
    for (int idx = 0;idx<mydat.nelem(); ++idx){
        // for (int idx2=0; idx2<mydat.nelem(); ++idx2){
        //     Mlarge(idx, idx2 + mydat.nelem()) = 0.0;
        // }
        Mlarge(idx, idx + mydat.nelem()) = 1.0;
    }
    a2corr = a2;
    for (int idx = 0; idx<mydat.nelem(); ++idx){
        a2corr(idx,idx) = a2(idx,idx) + 1.0;
    }
    a2inv = a2corr.inverse();
    Mlarge.block(mydat.nelem(),0,mydat.nelem(), mydat.nelem()) = - a2inv * a0;
    Mlarge.block(mydat.nelem(),mydat.nelem(),mydat.nelem(), mydat.nelem()) = - a2inv * a1;

    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    std::cout << "Inverse time: "
              << duration.count() / 1000000.0 << " seconds" << std::endl;

    ////////////////////////////////////////////////////////////////////////// 
    // Form RHS
    fRHS.resize(2 * mydat.nelem());
    fRHS.block(0,0,mydat.nelem(), 1) = Eigen::VectorXcd::Zero(mydat.nelem());
    fRHS.block(mydat.nelem(), 0, mydat.nelem(), 1) = a2inv * mydat.vs();
    
    ////////////////////////////////////////////////////////////////////////// 
    //eigensolve

    start = high_resolution_clock::now();
    Eigen::ComplexEigenSolver<Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> > msolve;
    
    //  Eigen::EigenSolver<Eigen::MatrixXcd> eigensolver;
    msolve.compute(Mlarge);
    // std::cout << Mlarge.rows() << " " << Mlarge.cols() << std::endl;
    Eigen::Vector< std::complex<double>, Eigen::Dynamic>  eig_values = msolve.eigenvalues();
    Eigen::Matrix< std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>  eig_vecs = msolve.eigenvectors();
    
    std::cout << eig_values(0) << " " << eig_values(1) << std::endl;
    std::cout << eig_vecs.block(0,0,5,2) << std::endl;
    
    
    
    
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    std::cout << "Eigensolution time: "
              << duration.count() / 1000000.0 << " seconds" << std::endl;

    ////////////////////////////////////////////////////////////////////////// 
    //solution
    start = high_resolution_clock::now();
    vecr.resize(2*mydat.nelem(), mydat.vr().cols());
    vecr.block(0,0,mydat.nelem(), mydat.vr().cols()) = mydat.vr();
    vecr.block(mydat.nelem(),0,mydat.nelem(), mydat.vr().cols()) = Eigen::MatrixXcd::Zero(mydat.nelem(), mydat.vr().cols());
    vecref = vecr.transpose() * Mlarge * eig_vecs;
    eiginv = eig_vecs.inverse();
    vecf = eiginv * fRHS;
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    std::cout << "Rest of calculation time: "
              << duration.count() / 1000000.0 << " seconds" << std::endl;

    std::cout << "Hello" << std::endl;
    std::cout << "nrows, ncols: " << vecref.rows() << " " << vecref.cols() << std::endl;

    //evaluating at times
    start = high_resolution_clock::now();
    std::vector<double> vec_t = myfreq.t();
    Eigen::MatrixXcd mat_tout;
    
    mat_tout.resize(mydat.nelem(), myfreq.nt());
    Eigen::MatrixXcd mat_tmp;
    mat_tmp.resize(2 * mydat.nelem(),1);
    std::complex<double> multfact;
    for (int idx = 0; idx< myfreq.nt(); ++idx){
        for (int idx2 = 0;idx2<2*mydat.nelem();++idx2){
            if (eig_values(idx2).real() > 0){
                multfact = 0.0;
                multfact = std::exp(eig_values(idx2) * static_cast<std::complex<double> > (vec_t[idx]));
            } else {
                multfact = std::exp(eig_values(idx2) * static_cast<std::complex<double> > (vec_t[idx]));
            }
            mat_tmp(idx2,0)  = multfact * vecf(idx2);
        }
        // if (idx == 1){
        //     std::cout << mat_tmp;
        // }
        

        mat_tout.block(0,idx,mydat.nelem2(),1) = vecref * mat_tmp;
    }
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    std::cout << "Rest of calculation time: "
              << duration.count() / 1000000.0 << " seconds" << std::endl;

    // std::cout << "nt is: " << myfreq.nt() << std::endl;
    // std::cout << "nt*2*nelem*nelem2 is: " << myfreq.nt() *2 * mydat.nelem() * mydat.nelem2()<< std::endl;

    std::cout << mat_tout.block(0,0,mydat.nelem2(),5) << std::endl;
    
    //equiv calc of raw spectra
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> rawspec;
    Eigen::MatrixXd tseis;
    tseis = mat_tout.real();
    std::cout << tseis.block(0,0,mydat.nelem2(),5) << std::endl;
    std::cout << tseis.block(0,myfreq.nt()-4,mydat.nelem2(),4) << std::endl;

    // std::cout << std::exp(eig_values(0) * static_cast<std::complex<double> > (vec_t[1])) << std::endl;
    // std::cout << "size: " << eig_values.size() << std::endl;

    // for (int idx = 0; idx < myfreq.nt();++idx){
    //     if (std::isnan(vec_t[idx])){
    //         std::cout << "At idx: " << idx << std::endl;
    //         break;
    //     }
    // }
    // rawspec = modespectrafunctions::rawspectra(myfreq, mydat, soltol);
    

    // Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> outspec(
    //     mydat.nelem2(), myfreq.nt() / 2 + 1);

    rawspec = processfunctions::rawtime2freq(tseis,myfreq.nt(), myfreq.dt());
    std::cout << rawspec.block(0,0,mydat.nelem2(),2) << std::endl;
    // return outspec;

    // finding seismogram, can use same as for test:
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> raw_seis;
    raw_seis =
        processfunctions::filtfreq2time(
        rawspec, myfreq.df(), myfreq.f1(), myfreq.f2(), myfreq.nt(),
        0.0, myfreq.dt(), myfreq.tout());
    std::cout << "Hello\n";
    std::cout << "Hello, it is: " << myfreq.nt0() << std::endl;
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>
        fin_spec;
    fin_spec =
        modespectrafunctions::calc_fspectra(raw_seis, myfreq, mydat.nelem2());
// std::cout << "Hello again\n";


    // outputting to files
    std::ofstream myfile;
    std::string outputfilebase = "fspectra.fullcouple.r";
    std::string outputfilename;

    // spectra
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

    

    // seismogram
    outputfilebase = "tspectra.fullcouple.r";
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
    

    // }
return 0;
}