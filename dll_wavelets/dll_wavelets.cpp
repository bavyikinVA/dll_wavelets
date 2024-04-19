#include "dll_wavelets.h"
#include <cmath>
#include <iostream>

extern "C" {
    void morlet_wavelet(int num_scales, int rows, int cols, double* data, double* scales, double* result) {
        for (int scale = 0; scale < num_scales; ++scale) {
            for (int row = 0; row < rows; ++row) {
                for (int col = 0; col < cols; ++col) {
                    double w0 = 0;
                    for (int k = 0; k < cols; ++k) {
                        double t = (double)(k - col) / scales[scale];
                        w0 += data[row * cols + col] * 0.75 * std::exp(-(t * t) / 2) * std::cos(6 * t);
                    }
                    result[scale * rows * cols + row * cols + col] = w0 / std::sqrt(scales[scale]);
                }
            }
        }
    }
}

