#pragma once

#include <fstream>
#include "SolVector.hpp"
#include "../Common/Config.hpp"
#include "../DualGrid/Grid.hpp"
#include "BLAS.hpp"
#include "../Common/AD.hpp"

class SolMatrix {
    public:
        SolMatrix();
        SolMatrix(const Grid *grid);
        ~SolMatrix();

        void Init(const Grid *grid);
        void Zeroes(void);

        const zdouble* GetBlock(const unsigned long row, const unsigned long col) const;
        zdouble* GetBlock(const unsigned long row, const unsigned long col);
        void AddBlock(const unsigned long row, const unsigned long col, const zdouble* block, const zdouble alpha = 1);
        void SubtractBlock(const unsigned long row, const unsigned long col, const zdouble* block, const zdouble alpha = 1);
        void AddToDiag(const unsigned long row, const unsigned long col, const zdouble val);
        void MatrixVecMult(const SolVector &b, SolVector &x) const;
        void RowProduct(const SolVector &b, SolVector &x, const unsigned long row) const;

        const unsigned long* GetRowPtr(void) const;
        const unsigned long* GetColInd(void) const;

        void WriteValues(void);

    private:
        zdouble *val;
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
