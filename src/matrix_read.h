#ifndef MATRIX_READ_GUARD_H
#define MATRIX_READ_GUARD_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

class couplematrix {
   public:
    using vecint = Eigen::Matrix<int, Eigen::Dynamic, 1>;
    using veccd = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>;
    using matcd =
        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>;

    // constructor
    couplematrix(std::string, std::string, std::string);

    // access
    //     const vecint& ll()const { return m_ll; }

    // const veccd& ww() const{ return m_ww; }
    // const veccd& vs()const { return m_vs; }

    //     const matcd &a0() const { return m_a0; }
    //     const matcd&  a1() const { return m_a1; }
    //     const matcd& a2() const{ return m_a2; }
    //     const matcd& vr() const{ return m_vr; }

    //     const std::complex<double>& vs(int idx) const{ return m_vs(idx); }

    //     const int &nelem() const { return m_nelem; }
    //     const int& nelem2() const{ return m_nelem2; }
    //     const int& mtot() const{ return m_mtot; }
    vecint ll() const { return m_ll; }

    veccd ww() const { return m_ww; }
    veccd vs() const { return m_vs; }

    matcd a0() const { return m_a0; }
    matcd a1() const { return m_a1; }
    matcd a2() const { return m_a2; }
    matcd vr() const { return m_vr; }

    std::complex<double> vs(int idx) const { return m_vs(idx); }

    int nelem() const { return m_nelem; }
    int nelem2() const { return m_nelem2; }
    int mtot() const { return m_mtot; }

   private:
    vecint m_ll;
    veccd m_vs, m_ww;
    matcd m_a0, m_a1, m_a2, m_vr;
    int m_nelem, m_nelem2, m_mtot;
};

couplematrix::couplematrix(std::string filepath, std::string filePath2,
                           std::string filePath3) {
    //////////////////////////////////////////////////////////////////////////
    /////////////////////////    m_a0, m_a1, m_a2    /////////////////////////
    //////////////////////////////////////////////////////////////////////////

    // reading in
    std::ifstream infile(filepath, std::ifstream::binary);

    // check opened correctly
    if (!infile) {
        std::cout << "Cannot open matrix file!" << std::endl;
    }

    // find the length of the file
    infile.seekg(0, std::ios::end);
    // long int file_size = infile.tellg();
    std::size_t file_size = infile.tellg();
    infile.seekg(0, infile.beg);   // reset back to beginning

    // finding number of elements
    m_nelem = sqrt((file_size - 24) / (3 * 16));

    std::cout << file_size << std::endl;
    std::cout << m_nelem << std::endl;
    // resizing matrices using size found
    m_a0.resize(m_nelem, m_nelem);
    std::cout << m_nelem << std::endl;
    m_a1.resize(m_nelem, m_nelem);
    m_a2.resize(m_nelem, m_nelem);

    // placeholder byte size 4, ie same as integer
    int i;
    // long int matbytes;
    std::size_t matbytes;
    infile.read(reinterpret_cast<char *>(&i),
                sizeof(i));   // head placeholder
    double matbreak;          // placeholder between matrices

    // size of matrices in bytes
    // matbytes = static_cast<long int>(m_nelem) * static_cast<long
    // int>(m_nelem) * 16;
    matbytes = m_nelem * m_nelem * 16;
    std::cout << matbytes << std::endl;
    // read matrices
    infile.read(reinterpret_cast<char *>(m_a0.data()), matbytes);
    infile.read(reinterpret_cast<char *>(&matbreak), 8);   // placeholder
    infile.read(reinterpret_cast<char *>(m_a1.data()), matbytes);
    infile.read(reinterpret_cast<char *>(&matbreak), 8);   // placeholder
    infile.read(reinterpret_cast<char *>(m_a2.data()), matbytes);
    infile.read(reinterpret_cast<char *>(&i), sizeof(i));   // placeholder

    // close
    infile.close();
    // check
    if (!infile.good()) {
        std::cout << "Matrix file error occurred at reading time!" << std::endl;
    }

    //////////////////////////////////////////////////////////////////////////
    ////////////////////     reading in m_vs and m_vr     ////////////////////
    //////////////////////////////////////////////////////////////////////////

    std::ifstream vecfile(filePath2, std::ifstream::binary);
    // check opened correctly
    if (!vecfile) {
        std::cout << "Cannot open source and receiver vector file!"
                  << std::endl;
    }

    // find the length of the file
    vecfile.seekg(0, std::ios::end);
    std::size_t file_size2 = vecfile.tellg();
    vecfile.seekg(0, vecfile.beg);   // reset back to beginning

    // finding number of elements
    m_nelem2 = (static_cast<int>(file_size2) - 16) / (16 * m_nelem) - 1;
    std::cout << file_size2 << std::endl;
    std::cout << m_nelem2 << std::endl;
    // resizing matrices using size found
    m_vs.resize(m_nelem);
    m_vr.resize(m_nelem, m_nelem2);

    // reading
    vecfile.read(reinterpret_cast<char *>(&i), sizeof(i));   // head

    // size
    int vecbytes = 16 * m_nelem;

    // read matrices
    vecfile.read(reinterpret_cast<char *>(m_vs.data()), vecbytes);
    vecfile.read(reinterpret_cast<char *>(&i), sizeof(i));   // placeholder
    vecfile.read(reinterpret_cast<char *>(&i), sizeof(i));   // placeholder
    vecfile.read(reinterpret_cast<char *>(m_vr.data()), vecbytes * m_nelem2);

    // close
    vecfile.close();
    // check
    if (!vecfile.good()) {
        std::cout
            << "Source and receiver vector file error occurred at reading time!"
            << std::endl;
    }

    //////////////////////////////////////////////////////////////////////////
    ////////////////////    reading in m_ll and m_ww      ////////////////////
    //////////////////////////////////////////////////////////////////////////
    std::ifstream freqfile(filePath3, std::ifstream::binary);

    // check opened correctly
    if (!freqfile) {
        std::cout << "Cannot open mode file!" << std::endl;
    }

    // find the length of the file
    freqfile.seekg(0, std::ios::end);
    std::size_t file_size3 = freqfile.tellg();
    freqfile.seekg(0, freqfile.beg);   // reset back to beginning

    // finding number of elements
    m_mtot = (file_size3 - 16) / (20);
    // resizing matrices using size found
    m_ll.resize(m_mtot);
    m_ww.resize(m_mtot);

    // placeholder
    freqfile.read(reinterpret_cast<char *>(&i), sizeof(i));   // head
    std::cout << "Hello\n";
    // read matrices
    freqfile.read(reinterpret_cast<char *>(m_ww.data()), 16 * m_mtot);
    freqfile.read(reinterpret_cast<char *>(&i), sizeof(i));   // placeholder
    freqfile.read(reinterpret_cast<char *>(&i), sizeof(i));   // placeholder
    freqfile.read(reinterpret_cast<char *>(m_ll.data()), 4 * m_mtot);

    // close
    freqfile.close();
    // check
    if (!freqfile.good()) {
        std::cout << "Mode file error occurred at reading time!" << std::endl;
    }
};
#endif