#pragma once

#include "SolVector.hpp"
#include "../Common/Config.hpp"
#include "BLAS.hpp"
#include "../Common/AD.hpp"

class DiagonalMatrix{
    public:
        DiagonalMatrix();
        ~DiagonalMatrix();

        void Zeroes(void);
        zdouble* GetBlock(const unsigned long iBlock);
        const zdouble* GetBlock(const unsigned long iBlock) const;
        void SetBlock(const zdouble *block, const unsigned long iBlock);
        void MatrixVecMult(const SolVector &vec, SolVector &prod) const;
    private:
        zdouble *val;
        unsigned long len, nBlk, blk_len, nVar;
};
