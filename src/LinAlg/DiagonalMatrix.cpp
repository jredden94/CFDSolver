#include "DiagonalMatrix.hpp"
#include "BLAS.hpp"

DiagonalMatrix::DiagonalMatrix() { 
    Config *config = &Config::GetConfig();
    nBlk = config->GetNumEqn();
    nVar = config->GetNumVars();
    blk_len = nVar * nVar;
    len = nBlk * blk_len;
    val = new zdouble[len];
}
DiagonalMatrix::~DiagonalMatrix() { delete[] val; }

void DiagonalMatrix::Zeroes(void) { 
    for (auto i = 0ul; i < len; i++) val[i] = 0;
}
zdouble* DiagonalMatrix::GetBlock(const unsigned long iBlock) { 
    return &val[iBlock * blk_len];
}
const zdouble* DiagonalMatrix::GetBlock(const unsigned long iBlock) const { 
    return &val[iBlock * blk_len];
}
void DiagonalMatrix::SetBlock(const zdouble *block, const unsigned long iBlock) { 
    zdouble *this_block = GetBlock(iBlock);
    for (auto i = 0ul; i < blk_len; i++) this_block[i] = block[i];
}
void DiagonalMatrix::MatrixVecMult(const SolVector &vec, SolVector &prod) const { 
    for (auto i = 0ul; i < nBlk; i++) {
        const zdouble *mat_blk = GetBlock(i);
        const zdouble *vec_blk = vec.GetBlock(i);
        zdouble *prod_blk = prod.GetBlock(i);
        blas::gemv(mat_blk, vec_blk, prod_blk, nVar, nVar);
    }
}
