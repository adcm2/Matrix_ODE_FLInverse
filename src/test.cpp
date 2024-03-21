#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif
// #ifndef EIGEN_USE_BLAS
// #define EIGEN_USE_BLAS
// #endif
// #ifndef EIGEN_USE_LAPACKE_STRICT
// #define EIGEN_USE_LAPACKE_STRICT
// #endif

// #include <math.h>
// #include <omp.h>
// #include <stdio.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <vector>

#include "BlockPreconditioner.h"
#include "Timer_Class.h"
#include "filter_header.h"
#include "freq_setup.h"
#include "matrix_read.h"
#include "matrix_replacement.h"
#include "spectra_central.h"

using namespace std::chrono;
int
main() {
    //////////////////////////////////////////////////////////////////////////
    // file path information
    std::string pathstring, filePath, filePath2, filePath3;
    pathstring = std::filesystem::current_path();
    filePath = pathstring + "/matrix.bin";
    filePath2 = pathstring + "/vector_sr.bin";
    filePath3 = pathstring + "/freq_sph.bin";
    //////////////////////////////////////////////////////////////////////////
    double f1, f2, dt, tout, df0, wtb, t1, t2, soltol;
    int qex;

    // inputting data for parameters of spectra
    std::cin >> f1 >> f2 >> dt >> tout >> df0 >> wtb >> t1 >> t2 >> soltol >>
        qex;
    // timer1.start();
    // getting setup of frequencies etc used in idsm
    freq_setup myfreq(f1, f2, dt, tout, df0, wtb, t1, t2, qex);
    // timer1.stop("Time taken to set up frequencies");
    //////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////
    // extracting coupling matrices from binary files
    Timer timer1;
    // timer1.start();

    // actual extraction
    couplematrix mydat(filePath, filePath2, filePath3);
    using MATRIX =
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>;
    using VECTOR = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>;
    // timer1.stop("Time taken to read in matrices");
    // checking matrix replacement
    // {
    //     std::complex<double> myi(0.0, 1.0);
    //     std::complex<double> winp = 0.01 - myi * 0.01;
    //     MATRIX mat_a0 = mydat.a0();
    //     MATRIX mat_a1 = mydat.a1();
    //     MATRIX mat_a2 = mydat.a2();
    //     MATRIX *_a0, *_a1, *_a2;
    //     _a0 = &mat_a0;
    //     _a1 = &mat_a1;
    //     _a2 = &mat_a2;
    //     timer1.start();
    //     MATRIX A(mydat.nelem(), mydat.nelem()), B, C;
    //     A = mat_a0 + winp * mat_a1 + winp * winp * mat_a2;

    //     // include diagonal component
    //     for (int idx = 0; idx < A.rows(); ++idx) {
    //         A(idx, idx) = A(idx, idx) - winp * winp;
    //     }

    //     MatrixReplaceFT<std::complex<double>> rep_A;
    //     rep_A.attachMatrices(mat_a0, mat_a1, mat_a2);
    //     rep_A.frequency(winp);
    //     Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> vrhs, x0,
    //         x1;   // declare solver, using diagonal preconditioner for moment
    //     vrhs.resize(mydat.nelem());
    //     x0.resize(mydat.nelem());
    //     x1.resize(mydat.nelem());

    //     // find rhs
    //     for (int idx = 0; idx < mydat.nelem(); ++idx) {
    //         vrhs(idx) = mydat.vs(idx) / (myi * winp);
    //     }

    //     Eigen::BiCGSTAB<
    //         Eigen::Matrix<std::complex<double>, Eigen::Dynamic,
    //         Eigen::Dynamic>,
    //         Eigen::BlockPreconditioner<std::complex<double>>>
    //         solver;

    //     // finding guess
    //     std::vector<int> vecidx;
    //     vecidx = randomfunctions::findindex(winp.real(), myfreq.wtb(),
    //                                         mydat.ll(), mydat.ww().real());

    //     solver.setTolerance(0.0001);
    //     solver.compute(A);
    //     solver.preconditioner().setblock(A, vecidx[0], vecidx[1]);
    //     // std::cout << "Hello 1\n";
    //     x0 = solver.preconditioner().solve(vrhs);
    //     // std::cout << "Hello 1.5\n";

    //     Eigen::BiCGSTAB<MatrixReplaceFT<std::complex<double>>,
    //                     Eigen::BlockMatFreePreconditioner<std::complex<double>>>
    //         solver2;
    //     // std::cout << "Hello 2\n";
    //     solver2.compute(rep_A);
    //     // std::cout << "Hello 3\n";
    //     solver2.preconditioner().setblock(rep_A, vecidx[0], vecidx[1]);
    //     x1 = solver2.preconditioner().solve(vrhs);
    //     // std::cout << "Hello 5\n";
    //     double tmpdiff;
    //     for (int i = 0; i < mydat.nelem(); ++i) {
    //         tmpdiff += std::abs(x1(i) - x0(i));
    //     }
    //     double tmp2 = std::sqrt(std::abs((x1 - x0).dot(x1 - x0)));
    //     double tmp3 = std::sqrt(std::abs((x1).dot(x1)));
    //     // std::cout << lhs1 - lhs2 << std::endl;
    //     std::cout << "Total difference: " << tmpdiff << std::endl;
    //     std::cout << "Relative difference: " << tmp2 / tmp3 << std::endl;
    // }

    // evaluate raw spectra
    timer1.start();
    auto rawspec = modespectrafunctions::rawspectra(myfreq, mydat, soltol);
    timer1.stop("Forming A");
    timer1.start();
    auto rawspec2 =
        modespectrafunctions::rawspectralowmemory(myfreq, mydat, soltol);
    timer1.stop("Not forming A");

    // comparing answers:
    {
        double tmp;
        double tmppsum;
        for (int i = 0; i < rawspec.cols(); ++i) {
            tmp += std::abs(rawspec(0, i) - rawspec2(0, i));
            tmppsum += std::abs(rawspec(0, i));
        }
        MATRIX tmpdiff = rawspec - rawspec2;

        // print
        // std::cout << tmpdiff.block(0, 0, 3, 10) << std::endl;
        // std::cout << rawspec.block(0, 0, 3, 10) << std::endl;
        auto rawlen = rawspec.cols();
        MATRIX tmprow = tmpdiff.block(0, 0, 1, rawlen);
        MATRIX tmpmat = tmprow * tmprow.transpose();
        MATRIX tmprow2 = rawspec.block(0, 0, 1, rawlen);
        MATRIX tmpmat2 = tmprow2 * tmprow2.transpose();
        // std::cout << "Rows: " << tmp3.rows() << ", columns: " << tmp3.cols()
        // << std::endl;
        double tmp2 = std::sqrt(std::abs(tmpmat(0, 0)));
        double tmpsum = std::sqrt(std::abs(tmpmat2(0, 0)));
        // double tmp3 =
        //     std::sqrt(std::abs(rawspec.block(0, 0, 1, rawlen) *
        //                        rawspec.block(0, 0, 1, rawlen).transpose()));

        std::cout << "Difference summed: " << tmp << std::endl;
        std::cout << "Difference summed: " << tmp / tmppsum << std::endl;
        std::cout << "Relative difference: " << tmp2 / tmpsum << std::endl;
    }
    // std::cout << "Hello!\n";
    // timer1.stop("Time taken to get raw spectra");

    // timer1.start();
    // find seismogram
    auto raw_seis =
        modespectrafunctions::calc_seismogram(rawspec, myfreq, mydat.nelem2());
    // timer1.stop("Time taken to get raw seismogram");

    // timer1.start();
    // find filtered spectra
    // Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic,
    //   Eigen::RowMajor>
    // fin_spec;
    auto fin_spec =
        modespectrafunctions::calc_fspectra(raw_seis, myfreq, mydat.nelem2());
    // timer1.stop("Time taken to get filtered spectrum");

    //////////////////////////////////////////////////////////////////////////
    // timer1.start();
    // outputting to files
    std::ofstream myfile;
    std::string outputfilebase = "fspectra.r";
    std::string outputfilename;

    // spectra
    for (int oidx = 0; oidx < mydat.nelem2(); ++oidx) {
        // name of output
        outputfilename = outputfilebase + std::to_string(oidx + 1) + ".out" +
                         ".q" + std::to_string(qex);

        // opening and writing to file
        myfile.open(outputfilename, std::ios::trunc);
        for (int idx = myfreq.i12(); idx < myfreq.i22(); ++idx) {
            myfile << std::setprecision(7) << myfreq.f2(idx) * 1000 << ";"
                   << fin_spec(oidx, idx).real() << ";"
                   << fin_spec(oidx, idx).imag() << ";"
                   << std::abs(fin_spec(oidx, idx)) << std::endl;
        }
        myfile.close();
    }

    // seismogram
    outputfilebase = "tspectra.r";
    for (int oidx = 0; oidx < mydat.nelem2(); ++oidx) {
        // name
        outputfilename = outputfilebase + std::to_string(oidx + 1) + ".out";

        // opening and writing to file
        myfile.open(outputfilename, std::ios::trunc);
        for (int idx = 0; idx < myfreq.nt(); ++idx) {
            // double tval = idx * mymode.dt / 3600;
            if (myfreq.t(idx) < myfreq.tout()) {
                myfile << std::setprecision(7) << myfreq.t(idx) << ";"
                       << raw_seis(oidx, idx) << std::endl;
            }
        }
        myfile.close();
    }
    // timer1.stop("Time to output data");
    return 0;
}
