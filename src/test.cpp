#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include "blockindex.h"
int
main() {
    Eigen::Matrix<int, Eigen::Dynamic, 1> ll;
    Eigen::Matrix<double, Eigen::Dynamic, 1> ww;

    // fill out
    // ll.resize(10);
    // ww.resize(10);
    // ll(0) = 0;
    // ll(1) = 1;
    // ll(2) = 1;
    // ll(3) = 1;
    // ll(4) = 2;
    // ll(5) = 2;
    // ll(6) = 2;
    // ll(7) = 2;
    // ll(8) = 2;
    // ll(9) = 2;
    // ww(0) = 0.0;
    // ww(1) = 1.0;
    // ww(2) = 1.0;
    // ww(3) = 2.0;
    // ww(4) = 2.0;
    // ww(5) = 2.0;
    // ww(6) = 3.0;
    // ww(7) = 4.0;
    // ww(8) = 5.0;
    // ww(9) = 6.0;
    ll.resize(6);
    ww.resize(6);
    ll(0) = 2;
    ll(1) = 2;
    ll(2) = 1;
    ll(3) = 3;
    ll(4) = 3;
    ll(5) = 4;
    ww(0) = 0.001943;
    ww(1) = 0.002382;
    ww(2) = 0.002538;
    ww(3) = 0.002944;
    ww(4) = 0.003683;
    ww(5) = 0.004066;
    std::vector<int> validx;
    validx = randomfunctions::findindex(0.002, 0.0003, ll, ww);
    std::cout << validx[0] << " " << validx[1] << std::endl;
    return 0;
}
