#ifndef MATRIX_HEADER_GUARD_H
#define MATRIX_HEADER_GUARD_H

#include <math.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <chrono>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using std::cout;
using std::endl;
using std::ifstream;
using std::ios;
using namespace std::chrono;

class modematrix {
   public:
    // constructor
    modematrix(std::string);

    // operator
    //  Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> operator()(
    //      std::complex<double>,
    //      Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>);

    // solve for spectra
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> finv(
        std::complex<double>);

   private:
    using myvector = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>;
    using mymatrix =
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>;
    mymatrix a0, a1, a2, vr;
    myvector vs;
    int nelem, nelem2;
};

// define constructor
modematrix::modematrix(std::string filepath) {
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
    //  cout << a0(0, 0) << endl;
    //  cout << a0(1, 0) << endl;
    //  cout << a1(0, 0) << endl;
    //  cout << a1(nelem - 1, nelem - 1) << endl;
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

    // ifstream infile("./vector_sr.bin", ifstream::binary);
    ifstream vecfile("/home/adcm2/raidam/coupling/work/vector_sr.bin",
                     ifstream::binary);
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
};

// define operator
Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>
modematrix::finv(std::complex<double> w) {
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> A(
        nelem, nelem);
    A = a0 + w * a1 + w * w * a2;
    Eigen::BiCGSTAB<
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>,
        Eigen::DiagonalPreconditioner<std::complex<double> > >
        solver;
    solver.compute(A);
    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> x;
    x = solver.solve(vs);
    return x;
};

#endif