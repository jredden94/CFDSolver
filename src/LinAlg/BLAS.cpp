#include "BLAS.hpp"
#include <iostream>

void blas::gemv(const zdouble *mat, const zdouble *vec_x, zdouble *vec_y, 
        const unsigned long m, const unsigned long n, const zdouble alpha, const zdouble beta, const bool scale) {
    for (auto i = 0ul; i < n; i++) {
        vec_y[i] = scale ? vec_y[i] * beta : zdouble{0.0};
        for (auto j = 0ul; j < m; j++) {
            vec_y[i] += alpha * mat[i * m + j] * vec_x[j]; 
        }
    }
}

void blas::gemm(const zdouble* mat_a, const zdouble *mat_b, zdouble *mat_c, 
        const unsigned long m) {
    for (auto i = 0ul; i < m; i++) {
        for (auto j = 0ul; j < m; j++) {
            mat_c[i*m+j] = 0.0;
            for (auto k = 0ul; k < m; k++) 
                mat_c[i*m+j] += mat_a[i*m+k] * mat_b[k*m+j];
        }
    }
}

void blas::InverseMatrix(const zdouble *mat, zdouble *mat_inv, const unsigned short m) {
    // Copy mat to mat2
    zdouble mat2[m * m];
    for (size_t i = 0; i < m * m; i++) mat2[i] = mat[i];

#define M(I,J) mat_inv[I * m + J]
#define A(I,J) mat2[I * m + J]

    /* Init mat_inv as identity */
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            M(i,j) = i == j ? 1.0 : 0.0;
        }
    }

    /* To Upper Triangular */
    for (auto i = 1ul; i < m; i++) {
        for (auto j = 0ul; j < i; j++) {
            zdouble weight = A(i,j) / A(j,j);
            for (auto k = j; k < m; k++) A(i,k) -= weight * A(j,k);
            for (auto k = 0ul; k <= j; k++) M(i,k) -= weight * M(j,k);
        }
    }

    /* Back Substitution */
    for (auto i = m; i > 0ul;) {
        i--;
        for (auto j = i + 1; j < m; j++)
            for (auto k = 0ul; k < m; k++) M(i,k) -= A(i,j) * M(j,k);

        for (auto k = 0ul; k < m; k++) M(i,k) /= A(i,i);
    }

#undef A
#undef M
}
