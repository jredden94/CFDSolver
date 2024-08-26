#pragma once

#include <fstream>
#include "SolVector.hpp"
#include "../Common/Config.hpp"
#include "../DualGrid/Grid.hpp"
#include "BLAS.hpp"

class SolMatrix {
    public:
        SolMatrix();
        SolMatrix(const Grid *grid);
        ~SolMatrix();

        void Init(const Grid *grid);
        void Zeroes(void);

        const double* GetBlock(const unsigned long row, const unsigned long col) const;
        double* GetBlock(const unsigned long row, const unsigned long col);
        void AddBlock(const unsigned long row, const unsigned long col, const double* block, const double alpha = 1);
        void SubtractBlock(const unsigned long row, const unsigned long col, const double* block, const double alpha = 1);
        void AddToDiag(const unsigned long row, const unsigned long col, const double val);
        void MatrixVecMult(const SolVector &b, SolVector &x) const;
        void RowProduct(const SolVector &b, SolVector &x, const unsigned long row) const;

        const unsigned long* GetRowPtr(void) const;
        const unsigned long* GetColInd(void) const;

        void WriteValues(void);

    private:
        double *val;
        unsigned long *row_ptr;
        unsigned long *col_ind;
        unsigned long nnz;
        unsigned long val_len;
        unsigned short blk_len;

        Config *config;
        const Grid *grid;
        unsigned short nVar;
        unsigned long nEqn;

        void Init(void);
};
