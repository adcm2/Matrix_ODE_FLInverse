// #ifndef EIGEN_DONT_PARALLELIZE
// #define EIGEN_DONT_PARALLELIZE
// #endif
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

#include "filter_header.h"
#include "freq_setup.h"
#include "matrix_read.h"
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
    // extracting coupling matrices from binary files
    auto start = high_resolution_clock::now();   // time start

    // actual extraction
    couplematrix mydat(filePath, filePath2, filePath3);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    std::cout << "Time taken to read in matrices: "
              << duration.count() / 1000000.0 << " seconds" << std::endl;

    //////////////////////////////////////////////////////////////////////////
    double f1, f2, dt, tout, df0, wtb, t1, t2, soltol;
    int qex;

    // inputting data for parameters of spectra
    std::cin >> f1 >> f2 >> dt >> tout >> df0 >> wtb >> t1 >> t2 >> soltol >>
        qex;

    // getting setup of frequencies etc used in idsm
    freq_setup myfreq(f1, f2, dt, tout, df0, wtb, t1, t2, qex);

    //////////////////////////////////////////////////////////////////////////
    // evaluate raw spectra
    auto rawspec = modespectrafunctions::rawspectra(myfreq, mydat, soltol);

    // find seismogram
    auto raw_seis =
        modespectrafunctions::calc_seismogram(rawspec, myfreq, mydat.nelem2());

    // find filtered spectra
    auto fin_spec =
        modespectrafunctions::calc_fspectra(raw_seis, myfreq, mydat.nelem2());

    //////////////////////////////////////////////////////////////////////////
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
            myfile << std::setprecision(17) << myfreq.f2(idx) * 1000 << ";"
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
                myfile << std::setprecision(17) << myfreq.t(idx) << ";"
                       << raw_seis(oidx, idx) << std::endl;
            }
        }
        myfile.close();
    }

    return 0;
}
