#ifndef BLOCKINDEX_GUARD_H
#define BLOCKINDEX_GUARD_H

#include <cstdlib>
#include <iostream>
#include <vector>

namespace randomfunctions {
std::vector<int>
findindex(const double& w, const double& wtb,
          const Eigen::Matrix<int, Eigen::Dynamic, 1>& ll,
          const Eigen::Matrix<double, Eigen::Dynamic, 1>& ww) {
    int maxidx = ll.rows() - 1;
    std::vector<int> vectmp = {0, 0};
    double wtb_tmp;
    double w_tmp;
    bool exitcond = false;

    // check positivity of wtb
    if (wtb < 0.0) {
        wtb_tmp = 0.0;
    } else {
        wtb_tmp = wtb;
    }

    // check w
    //  if (w < 0.0){
    //      w_tmp = 0.0;
    //  } else if (w > ww(maxidx)) {
    if (w + wtb < ww(0)) {
        return {0, 0};
    } else {
        w_tmp = w;
        int idxinc = 0;
        int idxstart = 0;
        for (int idx1 = 0; idx1 < ll.rows(); ++idx1) {
            if (std::abs(w - ww(idx1)) < wtb) {
                idxstart = idx1;
                break;
            }
            idxinc += 2 * ll(idx1) + 1;
            if (idx1 == ll.rows() - 1) {
                exitcond = true;
            }
        }
        if (exitcond) {
            vectmp[0] = 0;
            vectmp[1] = 0;
        } else {
            int idxinc2 = idxinc;
            for (int idx1 = idxstart; idx1 < ll.rows(); ++idx1) {
                if (std::abs(w - ww(idx1)) > wtb) {
                    // idxstart = idx1;
                    break;
                }
                idxinc2 += 2 * ll(idx1) + 1;
            }
            vectmp[0] = idxinc;
            vectmp[1] = idxinc2 - 1;
        }
        return vectmp;
    }
}

}   // namespace randomfunctions
// } else {
//     w_tmp = w;
// }

// check whether wtb > w
//  if (w_tmp - wtb_tmp <= 0){
//      vectmp[0] = 0;
//  } else {
//      auto upper = std::upper_bound(ww.begin(), ww.end(), w_tmp-wtb_tmp,
//      [](const double& x, const double& y){ return x <= y;}); int idx =
//      std::distance(ww.begin(), upper);
//      // std::cout << idx << std::endl;
//      vectmp[0] = std::min(std::max(0,idx), maxidx);
//  }

// check for w + wtb
//  if(w_tmp + wtb_tmp >= ww(maxidx)){
//      vectmp[1] = maxidx;
//  } else{
//      auto upper = std::upper_bound(ww.begin(), ww.end(), w_tmp + wtb_tmp);
//      int idx = std::distance(ww.begin(), upper);
//      // std::cout << idx << std::endl;
//      vectmp[1] = std::max(std::min(idx - 1,maxidx),0);
//  }

#endif