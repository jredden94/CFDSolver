#include "SolVector.hpp"
#include <fstream>

SolVector::SolVector() { 
    Config &config = Config::GetConfig();
    nVar = config.GetNumVars();
    nBlock = config.GetNumEqn();
    len = nVar * nBlock;
    Init();
}
SolVector::SolVector(const SolVector& vec) { 
    nVar = vec.nVar;
    nBlock = vec.nBlock;
    len = vec.len;
    Init();
}
SolVector::~SolVector() { 
    delete[] val;
    val = nullptr;
}

void SolVector::operator=(const SolVector &vec) { 
    for (auto i = 0ul; i < len; i++) val[i] = vec.val[i];
}

zdouble* SolVector::operator[](const unsigned long &iBlock) { return &val[iBlock * nVar]; }
const zdouble* const SolVector::operator[](const unsigned long &iBlock) const { return &val[iBlock * nVar]; }
void SolVector::operator+=(const SolVector &vec) { AddVector(vec); }
void SolVector::operator-=(const SolVector &vec) { SubtractVector(vec); }
void SolVector::operator*=(const zdouble scalar) {
    for (auto i = 0ul; i < len; i++) val[i] *= scalar;
}
void SolVector::operator/=(const zdouble scalar) {
    zdouble inv_scalar = 1.0 / scalar;
    for (auto i = 0ul; i < len; i++) val[i] *= inv_scalar;
}
void SolVector::PlusEqualScalarMult(const SolVector &vec, const zdouble scalar) {
    for (auto i = 0ul; i < len; i++) val[i] += scalar * vec.val[i];
}
void SolVector::MinusEqualScalarMult(const SolVector &vec, const zdouble scalar) {
    for (auto i = 0ul; i < len; i++) val[i] -= scalar * vec.val[i];
}

void SolVector::Init(void) { val = new zdouble[len]; }

void SolVector::AddBlock(const unsigned long& iNode, const vector<zdouble>& block, const zdouble& alpha) {
    for (unsigned short i = 0; i < nVar; i++) val[iNode * nVar + i] += block[i] * alpha;
}

void SolVector::AddBlock(const unsigned long& iNode, const zdouble* block, const zdouble& alpha) {
    for (unsigned short i = 0; i < nVar; i++) val[iNode * nVar + i] += block[i] * alpha;
}

void SolVector::SetBlock(const unsigned long& iNode, const vector<zdouble>& block, const zdouble& alpha) {
    for (unsigned short i = 0; i < nVar; i++) val[iNode * nVar + i] = block[i] * alpha;
}
void SolVector::SetBlock(const unsigned long& iNode, const zdouble* block, const zdouble& alpha) {
    for (unsigned short i = 0; i < nVar; i++) val[iNode * nVar + i] = block[i] * alpha;
}

const zdouble * const SolVector::GetBlock(const unsigned long& iNode) const { return &val[iNode * nVar]; }
zdouble* SolVector::GetBlock(const unsigned long& iNode) { return &val[iNode * nVar]; }
void SolVector::CopyBlock(const unsigned long& iBlock, zdouble *ret_block) { 
    const zdouble *block = GetBlock(iBlock);
    for (unsigned short i = 0; i < nVar; i++) ret_block[i] = block[i];
}

unsigned long SolVector::GetBlockCount(void) const { return nBlock; }

void SolVector::SubtractBlock(const unsigned long& iNode, const vector<zdouble>& block, const zdouble& alpha) { 
    for (unsigned short i = 0; i < nVar; i++) val[iNode * nVar + i] -= block[i] * alpha;
}

void SolVector::SubtractBlock(const unsigned long& iNode, const zdouble *block, const zdouble& alpha) { 
    for (unsigned short i = 0; i < nVar; i++) val[iNode * nVar + i] -= block[i] * alpha;
}

void SolVector::SetZeroes(void) { 
    for (auto i = 0ul; i < len; i++) val[i] = 0;
}

zdouble SolVector::DotProd(const SolVector& vec) const { 
    zdouble dot = 0;
    for (auto i = 0ul; i < len; i++) dot += val[i] * vec.val[i];
    return dot;
}

void SolVector::InitValues(const vector<zdouble>& x0) {
    assert(nVar == x0.size());
    for (auto i = 0ul; i < nBlock; i++) SetBlock(i, x0);
}

vector<zdouble> SolVector::ResNorm(void) const {
    vector<zdouble> res(nVar, 0);
    for (auto i = 0ul; i < nBlock; i++) {
        for (unsigned short j = 0; j < nVar; j++) {
            res[j] += pow(val[i * nVar + j], 2);
        }
    }

    for (unsigned short i = 0; i < nVar; i++) {
        res[i] = res[i] /  nBlock;
    }

    return res;
}
void SolVector::AddVector(const SolVector &vec, SolVector &result, const zdouble scale) const {
    for (auto i = 0ul; i < len; i++) result.val[i] = val[i] + scale * vec.val[i];
}

void SolVector::AddVector(const SolVector &vec, const zdouble scale) const {
    for (auto i = 0ul; i < len; i++) val[i] += scale * vec.val[i];
}

void SolVector::SubtractVector(const SolVector &vec, SolVector &result, const zdouble scale) const {
    const zdouble *vec_val = vec.val;
    for (auto i = 0ul; i < len; i++) result.val[i] = val[i] - scale * vec_val[i];
}

void SolVector::SubtractVector(const SolVector &vec, const zdouble scale) {
    for (auto i = 0ul; i < len; i++) val[i] -= scale * vec.val[i];
}

zdouble SolVector::ComputeMagnitude() const {
    zdouble mag = 0;
    for (auto i = 0ul; i < len; i++) mag += val[i] * val[i];
    return sqrt(mag);
}

zdouble SolVector::SquaredNorm(void) const { return DotProd(*this); }

zdouble SolVector::Norm(void) const { return sqrt(SquaredNorm() ); }

void SolVector::WriteValues(string filename) const {
    ofstream f(filename);
    for (auto i = 0ul; i < len; i++) f << val[i] << "\n";
    f.close();
}
