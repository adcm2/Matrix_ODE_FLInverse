#ifndef SPECTRA_EVAL_GUARD_H
#define SPECTRA_EVAL_GUARD_H

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

#include "BlockPreconditioner.h"
#include "blockindex.h"
#include "filter_header.h"

class modespectra {
   public:
    // constructor
    modespectra(std::string, std::string, std::string, double, double, double,
                double, double, double, double, double);

    // frequency of particular point
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>
    rawspectra();
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>
    fspectra();

    // solve for spectra
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> finv(
        const std::complex<double>&);

    // postprocessing
    Eigen::Matrix<double, 1, Eigen::Dynamic> postprocess(
        const Eigen::Matrix<std::complex<double>, 1, Eigen::Dynamic>&);
    Eigen::Matrix<std::complex<double>, 1, Eigen::Dynamic> postprocessf(
        const Eigen::Matrix<double, 1, Eigen::Dynamic>&);

    // number of points
    int nt, nt0;
    int i1, i2, mtot;
    int i12, i22;
    double df, ep, df2, dt;
    std::vector<double> w;
    std::vector<double> t;
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> fspec;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> tseis;
    int nelem, nelem2;

   private:
    double f1, f2, tout, df0, wtb, t1, t2;

    using myvector = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>;
    using mymatrix =
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>;
    mymatrix a0, a1, a2, vr;
    myvector vs;
    myvector ww;
    Eigen::Matrix<int, Eigen::Dynamic, 1> ll;
};

modespectra::modespectra(std::string filepath, std::string filePath2,
                         std::string filePath3, double f1, double f2, double dt,
                         double tout, double df0, double wtb, double t1,
                         double t2)
    : f1{f1 / 1000.0},
      f2{f2 / 1000.0},
      dt{dt},
      df0{df0 / 1000.0},
      tout{tout * 3600.0},
      t1{t1 * 3600.0},
      t2{std::min(t2, tout) * 3600.0},
      wtb{wtb * 3.1415926535 / 500.0} {
    /////////////////////////////////////////////////////////////////////////////////////
    // Setup
    /////////////////////////////////////////////////////////////////////////////////////
    // Nyquist frequency
    double fn = 0.5 / dt;

    // check for f2
    if (fn < this->f2) {
        std::cout << "f2 is greater than the Nyquist frequency for the time "
                     "step. Behaviour may be unexpected"
                  << std::endl;
        this->f2 = fn;
    };

    // multipliers for epsilon in inverse FL and q (need to clarify this)
    int mex = 5;
    int qex = 4;
    ep = mex / this->tout;
    // David's version
    //  df = ep / (6.28318530718 * qex);

    // my version
    df = 1.0 / this->tout;             // find df
    nt = std::ceil(1.0 / (df * dt));   // find nt
    int ne = static_cast<int>(log(static_cast<double>(nt)) / log(2.0) + 3);
    nt = pow(2, ne);   // finish increase in density
    std::cout << "nt: " << this->nt << std::endl;
    df = 1.0 / (nt * dt);   // new df

    // finding new values of f1 and f2 and corresponding integer numbers
    i1 = std::max(static_cast<int>(std::floor(this->f1 / df)), 0);
    i2 = static_cast<int>(std::floor(this->f2 / df)) + 2;
    this->f1 = (i1 - 1) * df;
    this->f2 = (i2 - 1) * df;

    // fill out w, remembering w = 2pi f, and only needing to go up to nt/2+1 as
    // half length of time series
    for (int idx = 0; idx < nt / 2 + 1; ++idx) {
        w.push_back(2.0 * 3.1415926535 * df * static_cast<double>(idx));
    }

    // fill out t
    for (int idx = 0; idx < nt; ++idx) {
        t.push_back(dt * static_cast<double>(idx));
    }
    this->nt0 = floor(1.0 / this->df0 * this->dt);
    if (this->nt0 > this->nt) {
        int ne2 = log(static_cast<double>(this->nt0)) / log(2.0) + 1;
        this->nt0 = pow(2, ne2);
    } else {
        this->nt0 = this->nt;
    }

    df2 = 1.0 / (nt0 * dt);
    i12 = std::max(static_cast<int>(floor(this->f1 / df2)), 1);
    i22 = static_cast<int>(floor(this->f2 / df2)) + 1;

    // std::cout << "nt: " << nt << ", nt0: " << nt0 << std::endl;
    //////////////////////////////////////////////////////////////////////////
    ///////////////      file opening and initialisation       ///////////////
    //////////////////////////////////////////////////////////////////////////
    // reading in for a0, a1, a2:
    ifstream infile(filepath, ifstream::binary);

    // check opened correctly
    if (!infile) {
        cout << "Cannot open file!" << endl;
    }

    // find the length of the file
    infile.seekg(0, ios::end);
    int file_size = infile.tellg();
    // cout << "Size of the file is " << file_size << " bytes" << endl;
    infile.seekg(0, infile.beg);   // reset back to beginning

    // finding number of elements
    nelem = sqrt((file_size - 24) / (3 * 16));
    // std::cout << nelem << std::endl;

    // placeholder byte size 4, ie same as integer
    int i, matbytes;
    infile.read(reinterpret_cast<char *>(&i), sizeof(i));   // head placeholder
    double matbreak;   // placeholder between matrices

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
    // std::cout << a0.cwiseAbs().block(0,0,10,10)/a0.cwiseAbs()(0,0) << std::endl;
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

    ifstream vecfile(filePath2, ifstream::binary);
    // check opened correctly
    if (!vecfile) {
        cout << "Cannot open file!" << endl;
    }

    // find the length of the file
    vecfile.seekg(0, ios::end);
    int file_size2 = vecfile.tellg();
    // cout << "Size of the file is " << file_size2 << " bytes" << endl;
    vecfile.seekg(0, vecfile.beg);   // reset back to beginning

    // finding number of elements
    nelem2 = (file_size2 - 16) / (16 * nelem) - 1;
// std::cout << nelem2 << std::endl;
    //////////////////////////////////////////////////////////////////////////
    ///////////////   reading in receiver and source vectors   ///////////////
    //////////////////////////////////////////////////////////////////////////
    vecfile.read(reinterpret_cast<char *>(&i), sizeof(i));   // head

    // resizing matrices using size found
    vs.resize(nelem);
    vr.resize(nelem, nelem2);

    // size
    int vecbytes = 16 * nelem;

    // read matrices
    vecfile.read(reinterpret_cast<char *>(vs.data()), vecbytes);
    vecfile.read(reinterpret_cast<char *>(&i), sizeof(i));   // placeholder
    vecfile.read(reinterpret_cast<char *>(&i), sizeof(i));   // placeholder
    vecfile.read(reinterpret_cast<char *>(vr.data()), vecbytes * nelem2);

    //////////////////////////////////////////////////////////////////////////
    ///////////////   closing and checking closed correctly    ///////////////
    //////////////////////////////////////////////////////////////////////////
    // close
    vecfile.close();
    // check
    if (!vecfile.good()) {
        cout << "Error occurred at reading time!" << endl;
    }

    //////////////////////////////////////////////////////////////////////////
    ///////////////      file opening and initialisation       ///////////////
    //////////////////////////////////////////////////////////////////////////
    ifstream freqfile(filePath3, ifstream::binary);
    // ifstream freqfile("../freq_sph.bin",
    //   ifstream::binary);
    // check opened correctly
    if (!freqfile) {
        cout << "Cannot open file!" << endl;
    }

    // find the length of the file
    freqfile.seekg(0, ios::end);
    int file_size3 = freqfile.tellg();
    // cout << "Size of the file is " << file_size3 << " bytes" << endl;
    freqfile.seekg(0, freqfile.beg);   // reset back to beginning

    // finding number of elements
    mtot = (file_size3 - 16) / (20);
    // std::cout << mtot << std::endl;

    //////////////////////////////////////////////////////////////////////////
    ///////////////     reading in data into a0, a1 and a2 ///////////////
    //////////////////////////////////////////////////////////////////////////
    // placeholder
    freqfile.read(reinterpret_cast<char *>(&i), sizeof(i));   // head

    // resizing matrices using size found
    ll.resize(mtot);
    ww.resize(mtot);

    // read matrices
    freqfile.read(reinterpret_cast<char *>(ww.data()), 16 * mtot);
    freqfile.read(reinterpret_cast<char *>(&i), sizeof(i));   // placeholder
    freqfile.read(reinterpret_cast<char *>(&i), sizeof(i));   // placeholder
    freqfile.read(reinterpret_cast<char *>(ll.data()), 4 * mtot);
    // std::cout << "Hello: " << mtot << std::endl;
    // std::cout << ww(0) << std::endl;
    // std::cout << ww(1) << std::endl;
    // std::cout << ww(2) << std::endl;
    // std::cout << ll(0) << std::endl;
    // std::cout << ll(1) << std::endl;
    // std::cout << ll(2) << std::endl;
    // std::cout << this->wtb << std::endl;
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

// define inverse function
Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>
modespectra::finv(const std::complex<double>& winp) {

    //coupling matrix
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> A(
        nelem, nelem);

    // declare value of A
    A = a0 + winp * a1 + winp * winp * a2;

    // include diagonal component
    for (int idx = 0; idx < A.rows(); ++idx) {
        A(idx, idx) = A(idx, idx) - winp * winp;
    }
    std::complex<double> myi(0.0, 1.0);   // imaginary number i

    // declare solver, using diagonal preconditioner for moment
    Eigen::BiCGSTAB<
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>,
        Eigen::BlockPreconditioner<std::complex<double> > >
        solver;
    //sparse trial
    // Eigen::SparseMatrix<std::complex<double> > S;
    // S = A.sparseView();
    //          Eigen::BiCGSTAB<
    //                 Eigen::SparseMatrix<std::complex<double> >,
    //                 Eigen::IncompleteLUT<std::complex<double> > >
    //                 solver;
                // Eigen::BiCGSTAB<
                //     Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>,
                //     Eigen::DiagonalPreconditioner<std::complex<double> > >
                //     solver;
    // indices
    std::vector<int> vecidx;
    vecidx = randomfunctions::findindex(winp.real(), wtb, ll, ww.real());

    solver.setTolerance(1.0 * std::pow(10.0, -5));
    solver.compute(A);
    // solver.compute(S);
    solver.preconditioner().setblock(A, vecidx[0], vecidx[1]);
    // solver.preconditioner().setblock(A, 0, 0);
    // solver.preconditioner().setDroptol(0.01);

    // solution and rhs
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> x, vrhs, x0;
    vrhs.resize(nelem);
    x0.resize(nelem);

    // find rhs
    for (int idx = 0; idx < nelem; ++idx) {
        vrhs(idx) = vs(idx) / (myi * winp);
        // x0(idx) = vrhs(idx) / A(idx, idx);
    }
    Eigen::BlockPreconditioner<std::complex<double> > incond;
    incond.compute(A);
    incond.setblock(A, vecidx[0], vecidx[1]);
    x0 = incond.solve(vrhs);

    // solve and return
    // x = solver.solve(vrhs);
    x = solver.solveWithGuess(vrhs, x0);
    // if (solver.iterations() > 3) {
    //     std::cout << "Frequency: " << winp.real() << ". Index:" << vecidx[0]
    //               << ", " << vecidx[1] << std::endl;
    //     std::cout << "#iterations:     " << solver.iterations() << std::endl;
    // }
    return x;
};

Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>
modespectra::rawspectra() {
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> tmp =
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic,
                      Eigen::Dynamic>::Zero(nelem2 + 1, nt / 2 + 1);
    // tmp.resize(nelem2 + 1, nt);
    std::complex<double> myi(0.0, 1.0);
    double oneovertwopi = 1.0 / (2.0 * 3.1415926535);
    std::complex<double> imep = static_cast<std::complex<double> >(ep);

    // go through frequencies:
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> vlhs;
    // #pragma omp parallel for
    for (int idx = 0; idx < nt / 2 + 1; ++idx) {
        tmp(0, idx) = w[idx] * oneovertwopi;   // frequency
        // std::complex<double> winp;             // imaginary shifted frequency
        // winp = w[idx] - myi * imep;

        // run through all frequencies and if between f1 and f2 compute
        if (idx > i1 - 1 && idx < i2 + 1) {
            // calculate x in Ax = b
            vlhs = modespectra::finv(w[idx] - myi * imep);

            // find acceleration response using receiver vectors
            tmp.block(1, idx, nelem2, 1) = -(w[idx] - myi * imep) *
                                           (w[idx] - myi * imep) *
                                           this->vr.transpose() * vlhs;
        };
    };
    return tmp;
};

Eigen::Matrix<double, 1, Eigen::Dynamic>
modespectra::postprocess(const Eigen::Matrix<std::complex<double>, 1, Eigen::Dynamic>&
                             rawspec) {   // do Fourier transform

    using Float = double;
    using Complex = std::complex<Float>;
    using RealVector = FFTWpp::vector<Float>;
    using ComplexVector = FFTWpp::vector<Complex>;
    using namespace std::complex_literals;

    // set up for FT

    RealVector testFL(nt), checkFL(nt);
    ComplexVector outFL(nt / 2 + 1);

    // Form the plans.
    auto flag = FFTWpp::Measure | FFTWpp::Estimate;

    auto outview = FFTWpp::MakeDataView1D(outFL);
    auto outview2 = FFTWpp::MakeDataView1D(checkFL);
    auto backward_plan = FFTWpp::Plan(outview, outview2, flag);
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> tmpspec;
    Eigen::Matrix<double, 1, Eigen::Dynamic> vecout;

    tmpspec.resize(this->nt / 2 + 1, 1);

    // set up Hann filter
    double fac = 0.1;
    double f11 = this->f1;
    double f22 = this->f2;
    double f12 = f11 + fac * (f22 - f11);
    double f21 = f22 - fac * (f22 - f11);

    // fill out for FT
    for (int idx = 0; idx < nt / 2 + 1; ++idx) {
        // actual frequency
        double finp;
        finp = static_cast<double>(idx) * this->df;
        outFL[idx] =
            rawspec(0, idx) * filters::hann(&finp, &f11, &f12, &f21, &f22);
    }

    // execute FT
    backward_plan.Execute();

    // do corrections
    auto myit2 = checkFL.begin();

    // output
    vecout.resize(1, nt);
    for (int idx = 0; idx < nt; ++idx) {
        double tinp;
        tinp = static_cast<double>(idx) * dt;
        if (tinp < this->tout) {
            vecout(0, idx) = myit2[idx] * exp(this->ep * t[idx]) * df;
        }
    }
    return vecout;
};

Eigen::Matrix<std::complex<double>, 1, Eigen::Dynamic>
modespectra::postprocessf(const Eigen::Matrix<double, 1, Eigen::Dynamic>& rawspec) {
    using Float = double;
    using Complex = std::complex<Float>;
    using RealVector = FFTWpp::vector<Float>;
    using ComplexVector = FFTWpp::vector<Complex>;
    using namespace std::complex_literals;

    RealVector dtmp(nt0);
    ComplexVector wdense(nt0 / 2 + 1);
    auto flag2 = FFTWpp::Measure | FFTWpp::Estimate;
    auto inView = FFTWpp::MakeDataView1D(dtmp);
    auto outView = FFTWpp::MakeDataView1D(wdense);

    auto forward_plan = FFTWpp::Plan(inView, outView, flag2);

    auto myit4 = dtmp.begin();

    //   parameters for time-domain filter
    double tfac = 0.5;
    double t11 = this->t1;
    double t22 = this->t2;
    double t12 = t11 + tfac * (t22 - t11);
    double t21 = t22 - tfac * (t22 - t11);
    for (int idx = 0; idx < nt; ++idx) {
        double tinp;
        tinp = static_cast<double>(idx) * dt;
        myit4[idx] =
            rawspec(0, idx) * filters::hann(&tinp, &t11, &t12, &t21, &t22);
    }
    if (nt0 > nt) {
        for (int idx = nt; idx < nt0; ++idx) {
            myit4[idx] = 0.0;
        }
    }

    forward_plan.Execute();

    auto myit5 = wdense.begin();
    Eigen::Matrix<std::complex<double>, 1, Eigen::Dynamic> tmp;
    tmp.resize(1, nt0 / 2 + 1);
    for (int idx = 0; idx < nt0 / 2 + 1; ++idx) {
        tmp(0, idx) = myit5[idx] * dt;
    }

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

    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> rawspec;
    rawspec = modespectra::rawspectra();

    // raw output
    std::ofstream myfile;

    myfile.open("frawspec.out", std::ios::trunc);
    for (int idx = i1 - 1; idx < i2 + 1; ++idx) {
        myfile << df * idx * 1000.0 << ";" << rawspec(1, idx).real() << ";"
               << rawspec(1, idx).imag() << ";" << std::abs(rawspec(1, idx))
               << std::endl;
    }
    myfile.close();

    // find seismogram
    this->tseis.resize(this->nelem2, nt);
    for (int idx = 0; idx < 3; ++idx) {
        this->tseis.block(idx, 0, 1, this->nt) =
            modespectra::postprocess(rawspec.block(1 + idx, 0, 1, nt / 2 + 1));
    }

    // frequency spectra
    this->fspec.resize(this->nelem2, nt0 / 2 + 1);
    for (int idx = 0; idx < 3; ++idx) {
        this->fspec.block(idx, 0, 1, this->nt0 / 2 + 1) =
            modespectra::postprocessf(this->tseis.block(idx, 0, 1, this->nt));
    }

    return this->fspec;
    // return 0;
};

#endif