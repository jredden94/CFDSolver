#include "BLAS.hpp"
#include "GuassSiedel.hpp"

GSS::GSS() { 
    Config &config = Config::GetConfig();
    vec_len = config.GetNumVars();
    mat_len = vec_len * vec_len;
    nEqn = config.GetNumEqn();
    lu_lower = new double[mat_len];
    lu_upper = new double[mat_len];
    b = new double[vec_len];
    y = new double[vec_len];
    x = new double[vec_len];
}
GSS::~GSS() { 
    delete[] lu_lower;
    delete[] lu_upper;
    delete[] b;
    delete[] y;
    delete[] x;
}

void GSS::Solve(const SolMatrix &mat, const SolVector &b_vec, SolVector &x_vec, const unsigned short nSweeps) {
    //x_vec.SetZeroes();
    const unsigned long *row_ptr = mat.GetRowPtr();
    const unsigned long *col_ind = mat.GetColInd();
    const double *diag, *off_diag;

    for (unsigned short sweep = 0; sweep < nSweeps; sweep++) {
        // Loop through matrix rows
        for (auto i = 0ul; i < nEqn; i++) {
            // b = b_vec block i
            const double* b_blk = b_vec.GetBlock(i);
            for (auto k = 0ul; k < vec_len; k++) b[k] = b_blk[k];
            diag = mat.GetBlock(i, i);

            // Loop through blocks on row i
            for (unsigned long j = row_ptr[i]; j < row_ptr[i+1]; j++) {
                unsigned long col = col_ind[j];
                if (col == i) continue;

                // b -= off_diag * x_vec block col
                const double *x_blk = x_vec.GetBlock(col);
                off_diag = mat.GetBlock(i, col);

                blas::gemv(off_diag, b_blk, b, vec_len, vec_len, -1, 1, true);
            }

            // Solve diag * x = b
            LUSolve(diag, b, x);

            x_vec.SetBlock(i, x);
        }
    }
}

void GSS::LUDecomposition(const double *block) {
    for (size_t i = 0; i < mat_len; i++) {
        lu_lower[i] = 0;
        lu_upper[i] = 0;
    }

    for (size_t i = 0; i < vec_len; i++) {
        // Upper
        for (size_t k = i; k < vec_len; k++) {
            double sum = 0;
            for (size_t j = 0; j < i; j++) {
                sum += lu_lower[i * vec_len + j] * lu_upper[j * vec_len + k];
            }
            lu_upper[i * vec_len + k] = block[i * vec_len + k] - sum;
        }

        // Lower
        for (size_t k = i; k < vec_len; k++) {
            if (i == k) lu_lower[i * vec_len + i] = 1;
            else {
                double sum = 0;
                for (size_t j = 0; j < i; j++) {
                    sum += lu_lower[k * vec_len + j] * lu_upper[j * vec_len + i];
                }
                lu_lower[k * vec_len + i] = (block[k * vec_len + i] - sum) / lu_upper[i * vec_len + i];
            }
        }
    }
}

void GSS::LUSolve(const double *mat, const double *b, double *x) {
    LUDecomposition(mat);

    // Ly = b
    for (unsigned short i = 0; i < vec_len; i++) {
        y[i] = b[i];
        for (unsigned short j = 0; j < i; j++) {
            y[i] -= lu_lower[i*vec_len + j] * y[j];
        }
        y[i] /= lu_lower[i * vec_len + i];
    }

    // Ux = y
    for (int i = vec_len - 1; i >= 0; i--) {
        x[i] = y[i];
        for (unsigned short j = i + 1; j < vec_len; j++) {
            x[i] -= lu_upper[i*vec_len + j] * x[j];
        }
        x[i] /= lu_upper[i * vec_len + i];
    }

    /*
       cout << "\n-----------------------------------\n";
       cout << "\nMatrix\n";
    for (auto i = 0ul; i < vec_len; i++) {
        for (auto j = 0ul; j < vec_len; j++) {
            cout << mat[i * vec_len + j] << "\t";
        }
        cout << "\n";
    }

    cout << "\n\nx\n";
    for (auto i = 0ul; i < vec_len; i++) cout << x[i] << "\t";

    cout << "\n\nb\n";
    for (auto i = 0ul; i < vec_len; i++) cout << b[i] << "\t";
    cout << "\n-----------------------------------\n";
    */
}
