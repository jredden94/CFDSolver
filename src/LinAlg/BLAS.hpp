#pragma once

#include "../Common/AD.hpp"

namespace blas {
    void gemv(const zdouble *mat, const zdouble *vec_x, zdouble *vec_y, 
            const unsigned long m, const unsigned long n, const zdouble alpha = 1, const zdouble beta = 0, const bool scale = false);

    void gemm(const zdouble* mat_a, const zdouble *mat_b, zdouble *mat_c, 
            const unsigned long m);

    void InverseMatrix(const zdouble *mat, zdouble *mat_inv, const unsigned short m);
}
