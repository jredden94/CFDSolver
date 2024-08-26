#include "SolMatrix.hpp"

SolMatrix::SolMatrix() { }
SolMatrix::SolMatrix( const Grid *grid ) : grid(grid) { Init(); }
SolMatrix::~SolMatrix() { 
    delete[] val;
    delete[] col_ind;
    delete[] row_ptr;
}

void SolMatrix::Init(const Grid *grid) {
    this->grid = grid;
    Init();
}

void SolMatrix::Init() {
    Config &config = Config::GetConfig();
    nVar = config.GetNumVars();
    nEqn = config.GetNumEqn();

    const vector<Node> &nodes = grid->Nodes();
    nnz = nodes.size();

    for (const Node &node : nodes) nnz += node.Neighbors().size();

    blk_len = nVar * nVar;
    val_len = nnz * blk_len;
    val = new double[val_len];
    col_ind = new unsigned long[nnz];
    row_ptr = new unsigned long[nEqn + 1];

    row_ptr[0] = 0;
    unsigned long colCount = 0;
    for (auto i = 0ul; i < nEqn; i++) {
        const Node &n = nodes[i];
        const vector<unsigned long> nbrs = n.Neighbors();
        row_ptr[i+1] = row_ptr[i] + nbrs.size() + 1;
        col_ind[colCount++] = i;
        for (const auto &nbr : nbrs) col_ind[colCount++] = nbr;
    }

    Zeroes();
}

void SolMatrix::Zeroes() {  
    for (auto i = 0ul; i < val_len; i++) val[i] = 0;
}

const double* SolMatrix::GetBlock(const unsigned long row, const unsigned long col) const {
    for (unsigned long i = row_ptr[row]; i < row_ptr[row+1]; i++)
        if (col_ind[i] == col) return &val[i * blk_len];
    return nullptr;
}

double* SolMatrix::GetBlock(const unsigned long row, const unsigned long col) {
    for (unsigned long i = row_ptr[row]; i < row_ptr[row+1]; i++)
        if (col_ind[i] == col) return &val[i * blk_len];
    return nullptr;
}

void SolMatrix::AddBlock(const unsigned long row, const unsigned long col, const double* block, const double alpha) {
    double *b = GetBlock(row, col);
    for (auto i = 0ul; i < blk_len; i++) { 
        b[i] += alpha * block[i];
    }
}

void SolMatrix::SubtractBlock(const unsigned long row, const unsigned long col, const double* block, const double alpha) {
    double *b = GetBlock(row, col);
    for (auto i = 0ul; i < blk_len; i++) b[i] -= alpha * block[i];
}

void SolMatrix::AddToDiag(const unsigned long row, const unsigned long col, const double val) {
    double *block = GetBlock(row, col);
    for (unsigned short i = 0; i < nVar; i++) block[i * nVar + i] += val;
}

void SolMatrix::MatrixVecMult(const SolVector &b, SolVector &x) const {
    x.SetZeroes();
    for (auto i = 0ul; i < nEqn; i++) RowProduct(b, x, i);
}

void SolMatrix::RowProduct(const SolVector &b, SolVector &x, const unsigned long row) const {
    for (auto i = row_ptr[row]; i < row_ptr[row+1]; i++) {
        const unsigned long &col = col_ind[i];
        const double *mat_blk = this->GetBlock(row, col);
        const double *vec_blk = b.GetBlock(col);
        double *prod_blk = x.GetBlock(col);
        blas::gemv(mat_blk, vec_blk, prod_blk, nVar, nVar, 1, 1, true);
    }
}

const unsigned long* SolMatrix::GetRowPtr(void) const { return row_ptr; }

const unsigned long* SolMatrix::GetColInd(void) const { return col_ind; }
