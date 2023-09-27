#ifndef SPECTRA_HEADER_GUARD_H
#define SPECTRA_HEADER_GUARD_H

#include <math.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <FFTWpp/All>
// #include "./FFTWpp/All"
#include <chrono>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "filter_header.h"

class modespectra {
   public:
    // constructor
    modespectra(std::string, std::string, double, double, double, double,
                double, double, double, double);

    // operator
    // double operator()(double);

    // frequency of particular point
    // std::complex<double> freq_value(int);
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>
    rawspectra();
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>
    fspectra();

    // solve for spectra
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> finv(
        std::complex<double>);

    // number of points
    int nt;
    int i1, i2, mtot;
    double df, ep;
    std::vector<double> w;
    // std::vector<int> ll;
    // std::vector<std::complex<double> > ww;

   private:
    double f1, f2, dt, tout, df0, wtb, t1, t2;

    using myvector = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>;
    using mymatrix =
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>;
    mymatrix a0, a1, a2, vr;
    myvector vs;
    myvector ww;
    Eigen::Matrix<int, Eigen::Dynamic, 1> ll;
    int nelem, nelem2;
};

modespectra::modespectra(std::string filepath, std::string filePath2, double f1,
                         double f2, double dt, double tout, double df0,
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
    if (fn < this->f2) {
        std::cout << "f2 is greater than the Nyquist frequency for the time "
                     "step. Behaviour may be unexpected"
                  << std::endl;
        this->f2 = fn;
    };
    int mex = 5;
    int qex = 4;
    // std::cout << modespectra::f1 << " " << modespectra::f2 << std::endl;
    ep = mex / this->tout;

    df = ep / (6.28318530718 * qex);
    // std::cout << "df: " << df << std::endl;
    // std::cout << "ep: " << ep << std::endl;
    // std::cout << "dt: " << dt << std::endl;
    // std::cout << "t/dt: " << modespectra::tout / dt << std::endl;

    nt = std::ceil(1.0 / (df * dt));
    std::cout << "nt: " << nt << std::endl;
    int ne = static_cast<int>(log(static_cast<double>(nt)) / log(2.0) + 1);
    std::cout << "ne: " << ne << std::endl;
    nt = pow(2, ne);

    df = 1.0 / (nt * dt);
    std::cout << "df: " << df << std::endl;

    i1 = std::max(static_cast<int>(std::floor(this->f1 / df)), 0);
    this->f1 = (i1 - 1) * df;
    // std::cout << "f1: " << df + this->f1 << std::endl;
    std::cout << "nt: " << nt << std::endl;
    std::cout << "f1: " << f1 << std::endl;
    // std::cout << this->f1 << " " << this->f2 << std::endl;
    i2 = static_cast<int>(std::floor(this->f2 / df)) + 2;

    this->f2 = (i2 - 1) * df;
    // std::cout << this->f1 << " " << this->f2 << std::endl;

    // fill out w
    for (int idx = 0; idx < nt; ++idx) {
        w.push_back(df * static_cast<double>(idx));
    }

    // change dt
    this->dt = 3.1415926535 / ((this->nt - 1) * this->df);
    //////////////////////////////////////////////////////////////////////////
    ///////////////      file opening and initialisation       ///////////////
    //////////////////////////////////////////////////////////////////////////

    ifstream infile(filepath, ifstream::binary);

    // check opened correctly
    if (!infile) {
        cout << "Cannot open file!" << endl;
    }

    // find the length of the file
    infile.seekg(0, ios::end);
    int file_size = infile.tellg();
    cout << "Size of the file is " << file_size << " bytes" << endl;
    infile.seekg(0, infile.beg);   // reset back to beginning

    // finding number of elements
    nelem = sqrt((file_size - 24) / (3 * 16));
    cout << nelem << endl;

    //////////////////////////////////////////////////////////////////////////
    ///////////////     reading in data into a0, a1 and a2     ///////////////
    //////////////////////////////////////////////////////////////////////////
    // placeholder byte size 4, ie same as integer
    int i, matbytes;
    infile.read(reinterpret_cast<char *>(&i), sizeof(i));   // head placeholder
    double matbreak;   // placeholder between matrices

    // reading into Eigen Matrix complex double. Need to have dynamic at start
    // as only find size at compile time
    // Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> a0,
    // a1, a2;

    // resizing matrices using size found
    a0.resize(nelem, nelem);
    a1.resize(nelem, nelem);
    a2.resize(nelem, nelem);

    // size of matrices in bytes
    matbytes = nelem * nelem * 16;

    // read matrices
    infile.read(reinterpret_cast<char *>(a0.data()), matbytes);
    infile.read(reinterpret_cast<char *>(&matbreak), 8);   // placeholder
    infile.read(reinterpret_cast<char *>(a1.data()), matbytes);
    infile.read(reinterpret_cast<char *>(&matbreak), 8);   // placeholder
    infile.read(reinterpret_cast<char *>(a2.data()), matbytes);
    infile.read(reinterpret_cast<char *>(&i), sizeof(i));   // placeholder

    //////////////////////////////////////////////////////////////////////////
    ///////////////                   checks                 ///////////////
    //////////////////////////////////////////////////////////////////////////
    // cout << a0(0, 0) << endl;
    // cout << a0(1, 0) << endl;
    // cout << a0(2, 0) << endl;
    // // cout << a0(1, 1) << endl;
    // cout << a1(0, 0) << endl;
    // cout << a1(1, 0) << endl;
    // cout << a1(2, 0) << endl;
    // // cout << a1(1, 1) << endl;
    // cout << a2(0, 0) << endl;
    // cout << a2(1, 0) << endl;
    // cout << a2(2, 0) << endl;
    // cout << a2(1, 1) << endl;
    //  cout << a2(0, 0) << endl;
    //  cout << a2(1, 0) << endl;

    //////////////////////////////////////////////////////////////////////////
    ///////////////   closing and checking closed correctly    ///////////////
    //////////////////////////////////////////////////////////////////////////
    // close
    infile.close();
    // check
    if (!infile.good()) {
        cout << "Error occurred at reading time!" << endl;
    }

    //////////////////////////////////////////////////////////////////////////
    ///////////////      file opening and initialisation       ///////////////
    //////////////////////////////////////////////////////////////////////////

    // ifstream vecfile("./vector_sr.bin", ifstream::binary);
    ifstream vecfile(filePath2, ifstream::binary);
    // check opened correctly
    if (!vecfile) {
        cout << "Cannot open file!" << endl;
    }

    // find the length of the file
    vecfile.seekg(0, ios::end);
    int file_size2 = vecfile.tellg();
    cout << "Size of the file is " << file_size2 << " bytes" << endl;
    vecfile.seekg(0, vecfile.beg);   // reset back to beginning

    // finding number of elements
    nelem2 = (file_size2 - 16) / (16 * nelem) - 1;
    cout << nelem2 << endl;

    //////////////////////////////////////////////////////////////////////////
    ///////////////     reading in data into a0, a1 and a2     ///////////////
    //////////////////////////////////////////////////////////////////////////
    // // placeholder byte size 4, ie same as integer
    // int i, matbytes;
    vecfile.read(reinterpret_cast<char *>(&i), sizeof(i));   // head
    // resizing matrices using size found
    vs.resize(nelem);
    vr.resize(nelem, nelem2);

    int vecbytes;
    vecbytes = 16 * nelem;

    // read matrices
    vecfile.read(reinterpret_cast<char *>(vs.data()), vecbytes);
    vecfile.read(reinterpret_cast<char *>(&i), sizeof(i));   // placeholder
    vecfile.read(reinterpret_cast<char *>(&i), sizeof(i));   // placeholder
    vecfile.read(reinterpret_cast<char *>(vr.data()), vecbytes * nelem2);

    //////////////////////////////////////////////////////////////////////////
    ///////////////                   checks                 ///////////////
    //////////////////////////////////////////////////////////////////////////
    // cout << vs(0) << endl;
    // cout << vs(1) << endl;
    // cout << vr(0, 0) << endl;
    // cout << vr(1, 0) << endl;

    //////////////////////////////////////////////////////////////////////////
    ///////////////   closing and checking closed correctly    ///////////////
    //////////////////////////////////////////////////////////////////////////
    // close
    vecfile.close();
    // check
    if (!vecfile.good()) {
        cout << "Error occurred at reading time!" << endl;
    }

    ifstream freqfile("/home/adcm2/raidam/coupling/work/freq_sph.bin",
                      ifstream::binary);
    // check opened correctly
    if (!freqfile) {
        cout << "Cannot open file!" << endl;
    }

    // find the length of the file
    freqfile.seekg(0, ios::end);
    int file_size3 = freqfile.tellg();
    cout << "Size of the file is " << file_size3 << " bytes" << endl;
    freqfile.seekg(0, freqfile.beg);   // reset back to beginning

    // finding number of elements
    mtot = (file_size3 - 16) / (20);
    // cout << mtot << endl;

    // //////////////////////////////////////////////////////////////////////////
    // ///////////////     reading in data into a0, a1 and a2 ///////////////
    // //////////////////////////////////////////////////////////////////////////
    // // // placeholder byte size 4, ie same as integer
    // // int i, matbytes;
    freqfile.read(reinterpret_cast<char *>(&i), sizeof(i));   // head
    // resizing matrices using size found
    ll.resize(mtot);
    ww.resize(mtot);

    // int vecbytes;
    // vecbytes = 16 * nelem;

    // // read matrices
    freqfile.read(reinterpret_cast<char *>(ww.data()), 16 * mtot);
    freqfile.read(reinterpret_cast<char *>(&i), sizeof(i));   // placeholder
    freqfile.read(reinterpret_cast<char *>(&i), sizeof(i));   // placeholder
    freqfile.read(reinterpret_cast<char *>(ll.data()), 4 * mtot);

    //////////////////////////////////////////////////////////////////////////
    ///////////////                   checks                 ///////////////
    //////////////////////////////////////////////////////////////////////////
    // cout << ll(0) << endl;
    // cout << ll(1) << endl;

    // ww = ww / 6.28318530718;
    // cout << ww << endl;
    // cout << ww(1) << endl;
    // cout << ww(2) << endl;
    // cout << ww(3) << endl;

    //////////////////////////////////////////////////////////////////////////
    ///////////////   closing and checking closed correctly    ///////////////
    //////////////////////////////////////////////////////////////////////////
    // close
    freqfile.close();
    // check
    if (!freqfile.good()) {
        cout << "Error occurred at reading time!" << endl;
    }
};

// define operator
Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>
modespectra::finv(std::complex<double> winp) {
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> A(
        nelem, nelem);

    A = a0 + winp * a1 + winp * winp * a2;
    // B = A.transpose();
    for (int idx = 0; idx < A.rows(); ++idx) {
        A(idx, idx) = A(idx, idx) - winp * winp;
    }
    std::complex<double> myi(0.0, 1.0);
    Eigen::BiCGSTAB<
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>,
        Eigen::DiagonalPreconditioner<std::complex<double> > >
        solver;
    // Eigen::FullPivLU<
    //     Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> >
    //     lusolve;
    solver.setTolerance(1.0 * std::pow(10, -7));
    solver.compute(A);
    // lusolve.compute(A);
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> x, vrhs;
    vrhs.resize(nelem);
    // do im = 1,mtot
    //        md = 2*ll(im)+1
    //        if(.not.check_mode(w,ww(im),wtb)) then
    //           do i = 1,md
    //              cfac = 1.0_sp/(ww(im)**2-w**2)
    //              xs(nindex+i) = cfac*vs(nindex+i)/(ii*w)
    //           end do
    //        end if
    //        nindex = nindex+md
    //     end do
    int nindex = 0;
    for (int idx = 0; idx < mtot; ++idx) {
        int md = 2 * ll(idx) + 1;
        for (int idx2 = 0; idx2 < md; ++idx2) {
            // std::complex<double> cfac = 1.0 / (ww(idx) * ww(idx) - winp *
            // winp); vrhs(nindex + idx2) = cfac * vs(nindex + idx2) / (myi *
            // winp); vrhs(nindex + idx2) = cfac * vs(nindex + idx2);
            vrhs(nindex + idx2) = vs(nindex + idx2) / (myi * winp);
        }
        nindex = nindex + md;
    }
    // x = lusolve.solve(vrhs);
    x = solver.solve(vrhs);
    return x;
};

Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>
modespectra::rawspectra() {
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> tmp;
    tmp.resize(nelem2 + 1, nt);
    std::complex<double> myi(0.0, 1.0);

    // go through frequencies:
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> vlhs;
    for (int idx = 0; idx < nt; ++idx) {
        for (int idx2 = 0; idx2 < nelem2 + 1; ++idx2) {
            tmp(idx2, idx) = 0.0;
        };
        tmp(0, idx) = w[idx];
        std::complex<double> winp;
        winp = 6.28318530718 * w[idx] -
               myi * static_cast<std::complex<double> >(ep);
        if (idx > i1 - 1 && idx < i2 + 1) {
            vlhs = modespectra::finv(winp);
            // for (int idx2 = 0;idx2<nelem2;++idx2){
            // tmp(idx2,idx) =vr.adjoint() * vlhs;
            // }
            tmp.block(1, idx, nelem2, 1) =
                -winp * winp * this->vr.transpose() * vlhs;
        };
        // std::cout << idx << " " << w[idx] << std::endl;
    };
    return tmp;
};

Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>
modespectra::fspectra() {
    // do Fourier transform
    using Float = double;
    using Complex = std::complex<Float>;
    using RealVector = FFTWpp::vector<Float>;
    using ComplexVector = FFTWpp::vector<Complex>;
    using namespace std::complex_literals;

    // set up for FT
    int n = 2 * (nt - 1);
    RealVector testFL(n), checkFL(n);
    ComplexVector outFL(nt);

    // Form the plans.
    auto flag = FFTWpp::Measure | FFTWpp::Estimate;

    auto backward_plan = FFTWpp::MakePlan1D(outFL, checkFL, flag);
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> rawspec,
        tmpspec;
    rawspec = modespectra::rawspectra();
    tmpspec.resize(this->nelem2 + 1, this->nt);

    // set up Hann filter
    double fac = 0.1;
    double f11 = this->f1;
    double f22 = this->f2;
    double f12 = f11 + fac * (f22 - f11);
    double f21 = f22 - fac * (f22 - f11);

    // do Hann filter
    for (int idx = 0; idx < this->nt - 1; ++idx) {
        double finp;
        finp = idx * this->df;
        tmpspec(0, idx) = rawspec(0, idx);
        tmpspec.block(1, idx, this->nelem2, 1) =
            rawspec.block(1, idx, this->nelem2, 1) *
            filters::hann(&finp, &f11, &f12, &f21, &f22);
    }

    // fill out for FT
    for (int idx = 0; idx < nt; ++idx) {
        outFL[idx] = (tmpspec(1, idx));
    }

    // execute FT

    backward_plan.Execute();

    // do corrections
    auto myit2 = checkFL.begin();
    for (int idx = 0; idx < n; ++idx) {
        myit2[idx] = myit2[idx] *
                     exp(this->ep * dt * static_cast<double>(idx)) * df *
                     static_cast<double>(nt - 1) /
                     (2.0 * 3.1415926535 * static_cast<Float>(n));
    }

    // output
    std::ofstream myfile;
    myfile.open("tspectratest.out", std::ios::trunc);
    for (int idx = 0; idx < n; ++idx) {
        myfile << this->dt * static_cast<double>(idx) << ";" << myit2[idx]
               << std::endl;
    }
    myfile.close();

    return tmpspec;
    // return 0;
};

#endif